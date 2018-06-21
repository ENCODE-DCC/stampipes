params.help = false
params.threads = 1
params.chunk_size = 16000000

params.UMI = false
params.trim_to = 0
params.genome = ""
params.outdir = "output"

def helpMessage() {
  log.info"""
    Usage: nextflow run basic.nf \\
             --genome /path/to/genome \\
             --bams '1.bam,2.bam...' \\

  """.stripIndent();
}

dataDir = "$baseDir/../../../data"
genome_name = file(params.genome).baseName

bams = Channel.fromPath(params.bams).collect()

process merge {

  input:
  file 'in*.bam' from bams

  output:
  file 'merged.bam' into merged

  publishDir params.outdir

  script:
  """
  echo hi
  samtools merge 'merged.bam' in*.bam
  """
}

// TODO: single end
process dups {
  publishDir params.outdir
  input:
  file(merged)

  output:
  file 'marked.bam' into marked_bam
  file 'MarkDuplicates.picard'

  script:
  if (params.UMI)
    cmd = "UmiAwareMarkDuplicatesWithMateCigar"
    extra = "UMI_TAG_NAME=XD"
  if (!params.UMI)
    cmd = "MarkDuplicatesWithMateCigar"
    extra = ""
  """
  picard RevertOriginalBaseQualitiesAndAddMateCigar \
    "INPUT=${merged}" OUTPUT=cigar.bam \
    VALIDATION_STRINGENCY=SILENT RESTORE_ORIGINAL_QUALITIES=false SORT_ORDER=coordinate MAX_RECORDS_TO_EXAMINE=0

  picard "${cmd}" \
      INPUT=cigar.bam OUTPUT=marked.bam \
      "${extra}" \
      METRICS_FILE=MarkDuplicates.picard ASSUME_SORTED=true VALIDATION_STRINGENCY=SILENT \
      READ_NAME_REGEX='[a-zA-Z0-9]+:[0-9]+:[a-zA-Z0-9]+:[0-9]+:([0-9]+):([0-9]+):([0-9]+).*'
  """
}

marked_bam.into { bam_for_counts; bam_for_adapter_counts; bam_for_filter }
process filter {
  input:
  file bam from bam_for_filter

  output:
  file "filtered.bam" into filtered_bam

  script:
  """
  samtools view -b -F 512 marked.bam > filtered.bam
  """
}
filtered_bam.into { bam_for_hotspot2; bam_for_spot_score; bam_for_cutcounts; bam_for_density; bam_for_inserts; bam_for_nuclear }

process filter_nuclear {
  input:
  file bam from bam_for_nuclear
  file nuclear_chroms from file("${params.genome}.nuclear.txt")

  output:
  file 'nuclear.bam' into nuclear_bam

  script:
  """
  samtools index "${bam}"
  cat "${nuclear_chroms}" \
  | xargs samtools view -b "${bam}" \
  > nuclear.bam
  """
}


process hotspot2 {

  publishDir "${params.outdir}"
  container "fwip/hotspot2:latest"

  input:
  file(marked_bam) from bam_for_hotspot2
  file(mappable) from file(params.mappable)
  file(chrom_sizes) from file(params.chrom_sizes)
  file(centers) from file(params.centers)


  output:
  file('peaks/filtered*')

  script:
  """
  hotspot2.sh -F 0.5 -p varWidth_20_default \
    -M "${mappable}" \
    -c "${chrom_sizes}" \
    -C "${centers}" \
    "${marked_bam}" \
    'peaks'
  """

}

process spot_score {
  publishDir params.outdir

  input:
  file(bam) from bam_for_spot_score
  file(mappable) from file("${dataDir}/annotations/${genome_name}.K36.mappable_only.bed")
  file(chromInfo) from file("${dataDir}/annotations/${genome_name}.chromInfo.bed")

  output:
  file 'subsample.spot.out'

  script:
  """
  # random sample
	samtools view -h -F 12 -f 3 "$bam" \
		| awk '{if( ! index(\$3, "chrM") && \$3 != "chrC" && \$3 != "random"){print}}' \
		| samtools view -uS - \
		> paired.bam
	bash \$STAMPIPES/scripts/bam/random_sample.sh paired.bam subsample.bam 5000000

  # hotspot
  bash \$STAMPIPES/scripts/SPOT/runhotspot.bash \
    \$HOTSPOT_DIR \
    \$PWD \
    \$PWD/subsample.bam \
    "${genome_name}" \
    36 \
    DNaseI
  """
}

process bam_counts {
  publishDir params.outdir

  input:
  file(bam) from bam_for_counts

  output:
  file('tagcounts.txt') into bam_counts

  script:
  """
  python3 \$STAMPIPES/scripts/bwa/bamcounts.py \
    "$bam" \
    tagcounts.txt
  """
}

process count_adapters {
  publishDir params.outdir

  input:
  file(bam) from bam_for_adapter_counts

  output:
  file('adapter.counts.txt') into adapter_counts

  script:
  """
  bash "\$STAMPIPES/scripts/bam/count_adapters.sh" "${bam}" > adapter.counts.txt
  """
}

process preseq {
  publishDir params.outdir
  input:
  file nuclear_bam

  output:
  file 'preseq.txt'
  file 'preseq_targets.txt'
  file 'dups.hist'

  script:
  """
  python3 \$STAMPIPES/scripts/bam/mark_dups.py -i "${nuclear_bam}" -o /dev/null --hist dups.hist
  preseq lc_extrap -hist dups.hist -extrap 1.001e9 -s 1e6 -v > preseq.txt

  # write out preseq targets
  bash "\$STAMPIPES/scripts/utility/preseq_targets.sh" preseq.txt preseq_targets.txt
  """
}

process cutcounts {

  publishDir params.outdir

  input:
  file(fai) from file("${params.genome}.fai")
  file(filtered_bam) from bam_for_cutcounts

  output:
  file('fragments.starch')
  file('cutcounts.starch')
  file('cutcounts.bw')
  file('cutcounts.bed.bgz')
  file('cutcounts.bed.bgz.tbi')

  script:
  """
  bam2bed --do-not-sort \
  < "$filtered_bam" \
  | awk -v cutfile=cuts.bed -v fragmentfile=fragments.bed -f \$STAMPIPES/scripts/bwa/aggregate/basic/cutfragments.awk

  sort-bed fragments.bed | starch - > fragments.starch
  sort-bed cuts.bed | starch - > cuts.starch

  unstarch cuts.starch \
  | cut -f1-3 \
  | bedops -m - \
  | awk '{ for(i = \$2; i < \$3; i += 1) { print \$1"\t"i"\t"i + 1 }}' \
  > allbase.tmp

  unstarch cuts.starch \
  | bedmap --echo --count --delim "\t" allbase.tmp - \
  | awk '{print \$1"\t"\$2"\t"\$3"\tid-"NR"\t"\$4}' \
  > cutcounts.tmp

  starch cutcounts.tmp > cutcounts.starch

  awk 'BEGIN{OFS="\t"}{print \$1,\$2,\$3,\$5}' cutcounts.tmp \
  | grep -v chrM \
  > wig.tmp

  wigToBigWig -clip wig.tmp "${fai}" cutcounts.bw

  # tabix
  unstarch cutcounts.starch | bgzip > cutcounts.bed.bgz
  tabix -p bed cutcounts.bed.bgz
  """
}

process density {

  publishDir params.outdir

  input:
  file filtered_bam from bam_for_density
  file chrom_bucket from file(params.chrom_bucket)

  output:
  file 'density.starch'

  script:
  """
  mkfifo density.bed

  bam2bed -d \
  < "${filtered_bam}" \
  | cut -f1-6 \
  | awk '{ if( \$6=="+" ){ s=\$2; e=\$2+1 } else { s=\$3-1; e=\$3 } print \$1 "\t" s "\t" e "\tid\t" 1 }' \
  | sort-bed - \
  > density.bed \
  &

  unstarch "${chrom_bucket}" \
  | bedmap --faster --echo --count --delim "\t" - density.bed \
  | awk -v binI=20 -v win=75 \
        'BEGIN{ halfBin=binI/2; shiftFactor=win-halfBin } { print \$1 "\t" \$2 + shiftFactor "\t" \$3-shiftFactor "\tid\t" i \$4}' \
  | starch - \
  > density.starch
  """

}

process insert_sizes {

  publishDir params.outdir

  input:
  file nuclear_bam from bam_for_inserts
  file nuclear_chroms from file("${params.genome}.nuclear.txt")

  output:
  file 'CollectInsertSizeMetrics.picard*'

  script:
  """
  picard CollectInsertSizeMetrics \
    "INPUT=${nuclear_bam}" \
    OUTPUT=CollectInsertSizeMetrics.picard \
    HISTOGRAM_FILE=CollectInsertSizeMetrics.picard.pdf \
    VALIDATION_STRINGENCY=LENIENT \
    ASSUME_SORTED=true

  cat CollectInsertSizeMetrics.picard \
  | awk '/## HISTOGRAM/{x=1;next}x' \
  | sed 1d \
  > hist.txt

  python3 "\$STAMPIPES/scripts/utility/picard_inserts_process.py" hist.txt > CollectInsertSizeMetrics.picard.info
  """
}