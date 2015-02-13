###################
# This makefile takes one set of FastQ files and aligns them using BWA, using a trimmer.
# It then provides a sorted BAM file.  The TMPDIR is best passed in as an empty directory
# to prevent collisions.
###################
# These variables must be passed in or set for the makefile to work.  If the genome's
# FAI file is not at $(BWAINDEX).fai, then it must also be specified under FAI.
###################
# BWAINDEX=/path/to/genome
# FASTQ1_FILE=/path/to/R1.fq.gz
# FASTQ2_FILE=/path/to/R2.fq.gz
# OUTBAM=/path/to/out.sorted.bam
# TRIMSTATS=/path/to/trimstats.txt
# ADAPTERFILE=/path/to/adapters.txt
###################
# REQUIRED MODULES
###################
# module load java
# module load bwa/0.7.0
# module load samtools
# module load python
###################

SHELL = bash

# NSLOTS is the default sge environment variable letting the process know
# how many threads it has reserved
# We can use that for the default number of threads to use for the bwa
# aln process
NSLOTS ?= 1
THREADS ?= $(NSLOTS)

# Use the default PATH version of each of these programs unless overridden
SAMTOOLS ?= samtools
BWA ?= bwa

FAI ?= $(BWAINDEX).fai

READ_LENGTH ?= 36

# Ideally this will be set in the environment or on the command line to an actual temp directory
TMPDIR ?= $(shell pwd)
ADAPTER_REPORT ?= $(shell pwd)/trimming_adapters.txt

# Discard reads failing chastity filter, use 32 bases as a seed
# and set the mismatch level appropriate to the read length
# NOTE: THESE ARE NOT WORKING FOR CHANGING IN QMAKE SCRIPT
BWA_ALN_OPTIONS ?= -Y -l 32 -n 0.04

# Include up to 10 matches for each read and
# set the inset size to a generous 750
BWA_SAMPE_OPTIONS ?= -n 10 -a 750

ADAPTER1 ?= P5_single_no_index
ADAPTER2 ?= P7_truseq_default # Probably at least one of these should not have a default

# List the intermediate files that can be deleted when we finish up
.INTERMEDIATE : 

all : info adapters align

info : 
	@echo "------"
	@echo "PROGRAM VERSIONS"
	@echo "------"
	@echo "bwa:"
	@$(BWA) 2>&1 | grep "Version"
	@echo "samtools:"
	@$(SAMTOOLS) 2>&1 | grep "Version"
	@echo "trim-adapters-illumina:"
	@trim-adapters-illumina -v
	@echo "------"
	@echo "METADATA"
	@echo "------"
	@echo "SAMPLE_NAME: " $(SAMPLE_NAME)
	@echo "BWAINDEX: " $(BWAINDEX)
	@echo "FAI: " $(FAI)
	@echo "READ LENGTH: " $(READ_LENGTH)
	@echo "ADAPTER1: " $(ADAPTER1)
	@echo "ADAPTER2: " $(ADAPTER2)
	@echo "ADAPTERFILE: " $(ADAPTERFILE)
	@echo "------"
	@echo "PROGRAM OPTIONS"
	@echo "------"
	@echo "BWA_ALN_OPTIONS: " $(BWA_ALN_OPTIONS)
	@echo "BWA_SAMPE_OPTIONS: " $(BWA_SAMPE_OPTIONS)
	@echo "------"

align : $(OUTBAM)

adapters: $(ADAPTER_REPORT)

# We don't really need the file written multiple times, but this'll solve races
$(ADAPTER_REPORT) : $(ADAPTERFILE)
	flock -x $@ awk '$$1 == "$(ADAPTER1)" {a1=$$0} $$1 == "$(ADAPTER2)" {a2=$$0} END {print a1; print a2}' $^ > $@

# Copy the final sorted bam to its finished place
$(OUTBAM) : $(TMPDIR)/align.sorted.bam
	rsync $^ $@

# Create sorted BAM files
# Note: samtools expects the output name without a .bam suffix and will
# add it, so take it away to prevent .bam.bam
$(TMPDIR)/align.sorted.bam : $(TMPDIR)/align.unsorted.bam
	time $(SAMTOOLS) sort $^ $(basename $@) && echo made $(TMPDIR)/align.sorted.bam >&2

# Create unsorted raw BAM files
$(TMPDIR)/align.unsorted.bam : $(TMPDIR)/align.sam
	time $(SAMTOOLS) view -bS -t $(FAI) $^ \
		| python $(STAMPIPES)/scripts/bwa/fix_bam_pairing.py - $@ \
		&& echo made $(TMPDIR)/align.unsorted.bam >&2

# Create the SAM files from each pair of SAI and FASTQ files
$(TMPDIR)/align.sam : $(TMPDIR)/R1.sai $(TMPDIR)/R2.sai $(TMPDIR)/trimmed.R1.fastq.gz $(TMPDIR)/trimmed.R2.fastq.gz
	time $(BWA) sampe $(BWA_SAMPE_OPTIONS) $(BWAINDEX) $^ > $@ && echo made $(TMPDIR)/align.sam >&2

# Make SAI alignments for each FASTQ file
$(TMPDIR)/%.sai : $(TMPDIR)/trimmed.%.fastq.gz
	time $(BWA) aln -t $(THREADS) $(BWA_ALN_OPTIONS) $(BWAINDEX) $(TMPDIR)/trimmed.$*.fastq.gz > $(TMPDIR)/$*.sai && echo made $(TMPDIR)/$*.sai >&2

# Trim each FASTQ file pair
# Keep track of the output by passing it to a file, so we can know how many pairs were trimmed
$(TMPDIR)/trimmed.R1.fastq.gz $(TMPDIR)/trimmed.R2.fastq.gz : $(FASTQ1_FILE) $(FASTQ2_FILE)
	time trim-adapters-illumina \
		-f $(ADAPTERFILE) \
		-1 $(ADAPTER1) -2 $(ADAPTER2) \
		$(FASTQ1_FILE) $(FASTQ2_FILE) \
		$(TMPDIR)/trimmed.R1.fastq.gz $(TMPDIR)/trimmed.R2.fastq.gz \
		&> $(TRIMSTATS) \
		&& echo trimmed $(SAMPLE_NAME) >&2
