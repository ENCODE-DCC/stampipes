import json
import os
import sys
import argparse
import logging
import requests
import subprocess

log_format = "%(asctime)s - %(name)s - %(levelname)s - %(message)s"

STAMPIPES = os.getenv('STAMPIPES', '~/stampipes')

script_options = {
    "quiet": False,
    "debug": False,
    "base_api_url": None,
    "token": None,
    "align_ids": [],
    "flowcell": None,
    "tag": None,
    "outfile": os.path.join(os.getcwd(), "run.bash"),
    "sample_script_basename": "run.bash",
    "qsub_prefix": ".proc",
    "dry_run": False,
    "no_mask": False,
    "bases_mask": None,
}

def parser_setup():

    parser = argparse.ArgumentParser()

    parser.add_argument("-q", "--quiet", dest="quiet", action="store_true",
        help="Don't print info messages to standard out.")
    parser.add_argument("-d", "--debug", dest="debug", action="store_true",
        help="Print all debug messages to standard out.")

    parser.add_argument("-a", "--api", dest="base_api_url",
        help="The base API url, if not the default live LIMS.")
    parser.add_argument("-t", "--token", dest="token",
        help="Your authentication token.  Required.")

    parser.add_argument("--alignment", dest="align_ids", type=int, action="append",
        help="Run for this particular alignment.")
    parser.add_argument("--flowcell", dest="flowcell_label",
        help="Run for this particular flowcell label.")
    parser.add_argument("--tag", dest="tag",
        help="Run for alignments tagged here.")

    parser.add_argument("-o", "--outfile", dest="outfile",
        help="Append commands to run this alignment to this file.")
    parser.add_argument("-b", "--sample-script-basename", dest="sample_script_basename",
        help="Name of the script that goes after the sample name.")
    parser.add_argument("--qsub-prefix", dest="qsub_prefix",
        help="Name of the qsub prefix in the qsub job name.  Use a . in front to make it non-cluttery.")
    parser.add_argument("-n", "--dry-run", dest="dry_run", action="store_true",
        help="Take no action, only print messages.")
    parser.add_argument("--no-mask", dest="no_mask", action="store_true",
        help="Don't use any barcode mask.")
    parser.add_argument("--bases_mask", dest="bases_mask",
        help="Set a bases mask.")

    parser.set_defaults( **script_options )
    parser.set_defaults( quiet=False, debug=False )

    return parser


class ProcessSetUp(object):

    def __init__(self, args, api_url, token): 

        self.token = token
        self.api_url = api_url
        self.qsub_scriptname = args.sample_script_basename
        self.qsub_prefix = args.qsub_prefix
        self.outfile = args.outfile
        self.dry_run = args.dry_run
        self.no_mask = args.no_mask
        self.headers = {'Authorization': "Token %s" % self.token}

    def api_single_result(self, url_addition=None, url=None):

        if url_addition:
           url = "%s/%s" % (self.api_url, url_addition)

        request = requests.get(url, headers=self.headers)

        if request.ok:
            logging.debug(request.json())
            return request.json()
        else:
            logging.error("Could not get data from %s" % url)
            logging.error(request)
            return None

    def api_list_result(self, url_addition=None, url=None):

        more = True
        results = []

        if url_addition:
            url = "%s/%s" % (self.api_url, url_addition)

        while more:

            logging.debug("Fetching more results for query %s" % url)

            request = requests.get(url, headers=self.headers)

            if not request.ok:
                logging.error(request)
                return None
            more_results = request.json()
            results.extend(more_results["results"])
            if more_results["next"]:
                url = more_results["next"]
            else:
                more = False

        return results

    def get_align_process_info(self, alignment_id):

        process_info = self.api_single_result("flowcell_lane_alignment/%d/processing_information" % alignment_id)

        if not process_info:
            logging.critical("Could not find processing info for alignment %d\n" % alignment_id)
            logging.critical(info)
            sys.exit(1)

        return process_info

    def get_process_template(self, align_id, process_template_id):

        if not process_template_id:
            logging.critical("No process template for alignment %d\n" % align_id)
            return None

        info = self.api_single_result("process_template/%d" % (process_template_id))

        if not info:
            logging.critical("Could not find processing template for ID %d\n" % process_template_id)
            sys.exit(1)

        return info

    def setup_alignment(self, align_id):

        processing_info = self.get_align_process_info(align_id)

        self.create_script(processing_info)

    def setup_tag(self, tag_slug):

        align_tags = self.api_list_result("tagged_object?content_type=47&tag__slug=%s" % tag_slug)

        [self.setup_alignment(align_tag["object_id"]) for align_tag in align_tags]

    def setup_flowcell(self, flowcell_label):

        logging.info("Setting up flowcell for %s" % flowcell_label)
        alignments = self.api_list_result("flowcell_lane_alignment?lane__flowcell__label=%s" % flowcell_label)

        [self.setup_alignment(alignment["id"]) for alignment in alignments]

    def add_script(self, align_id, processing_info, script_file, sample_name):

        if not self.outfile:
            logging.debug("Writing script to stdout")
            outfile = sys.stdout
        else:
            logging.debug("Logging script to %s" % self.outfile)
            outfile = open(self.outfile, 'a')

        outfile.write("cd %s && " % os.path.dirname(script_file))
        outfile.write("qsub -N %s%s-%s-ALIGN#%d -cwd -V -S /bin/bash %s\n\n" % (self.qsub_prefix, sample_name, processing_info['flowcell']['label'], align_id, script_file))

        outfile.close()

    def get_script_template(self, process_template):

        script_path = os.path.expandvars(process_template["process_version"]["script_location"])
        return open(script_path, 'r').read()

    def create_script(self, processing_info):

        lane = processing_info["libraries"][0]
        alignment = lane["alignments"][0]
        align_id = alignment["id"]

        if not "process_template" in alignment:
            logging.critical("Alignment %d has no process template" % align_id)
            return False

        process_template = self.get_process_template(align_id, alignment["process_template"])

        if not process_template:
            return

        flowcell_directory = processing_info['flowcell']['directory']

        if not flowcell_directory:
            logging.critical("Alignment %d has no flowcell directory for flowcell %s" % (align_id, processing_info['flowcell']['label']))
            return False

        fastq_directory = os.path.join(processing_info['flowcell']['directory'], "Project_%s" % lane['project'], "Sample_%s" % lane['samplesheet_name'])

        # Reset the alignment's sample name if we decied not to use the barcode index mask
        if self.no_mask:
            alignment['sample_name'] = "%s_%s_L00%d" % (lane['samplesheet_name'], lane['barcode_index'], lane['lane'])

        align_dir = "align_%d_%s_%s" % (alignment['id'], alignment['genome_index'], alignment['aligner'])
        if alignment['aligner_version']:
            align_dir = "%s-%s" % (align_dir, alignment['aligner_version'])

        script_directory = "%s/%s" % (fastq_directory, align_dir)

        script_file = os.path.join( script_directory, "%s-%s" % (alignment['sample_name'], self.qsub_scriptname) )
        logging.info(script_file)

        if self.dry_run:
            logging.info("Dry run, would have created: %s" % script_file)
            return

        if not os.path.exists(script_directory):
            logging.info("Creating directory %s" % script_directory)
            os.makedirs(script_directory)

        self.add_script(align_id, processing_info, script_file, alignment['sample_name'])

        outfile = open(script_file, 'w')
        outfile.write("set -e -o pipefail\n")
        outfile.write("export SAMPLE_NAME=%s\n" % alignment['sample_name'])
        outfile.write("export BWAINDEX=%s\n" % alignment['genome_index_location'])
        outfile.write("export GENOME=%s\n" % alignment['genome_index'])
        outfile.write("export ASSAY=%s\n" % lane['assay'])
        outfile.write("export READLENGTH=%s\n" % processing_info['flowcell']['read_length'])
        if processing_info['flowcell']['paired_end']:
            outfile.write("export PAIRED=True\n")
        outfile.write("export FLOWCELL_LANE_ID=%s\n" % lane['id'])
        outfile.write("export ALIGNMENT_ID=%s\n" % alignment['id'])
        outfile.write("export ALIGN_DIR=%s/%s\n" % (fastq_directory, align_dir))
        outfile.write("export R1_FASTQ=%s/%s_%s_R1.fastq.gz\n" %  (fastq_directory, processing_info['flowcell']['label'], alignment['sample_name']))
        outfile.write("export R2_FASTQ=%s/%s_%s_R2.fastq.gz\n" %  (fastq_directory, processing_info['flowcell']['label'], alignment['sample_name']))
        outfile.write("export FASTQ_DIR=%s\n" % fastq_directory)
        outfile.write("export FLOWCELL=%s\n" % processing_info['flowcell']['label'])
        if "barcode1" in lane and lane["barcode1"]:
            p7_adapter = lane['barcode1']['adapter7']
            p5_adapter = lane['barcode1']['adapter5']
            if "barcode2" in lane and lane['barcode2']:
                # Override the "default" end adapter from barcode1
                # TODO: Make sure we want adapter7, double-check lims methods
                p5_adapter = lane['barcode2']['adapter7']
            outfile.write("export ADAPTER_P7=%s\n" % p7_adapter)
            outfile.write("export ADAPTER_P5=%s\n" % p5_adapter)

            # Process with UMI if the barcode has one and this is a dual index
            # flowcell
            if lane['barcode1']['umi'] and processing_info['flowcell']['dual_index']:
                outfile.write("export UMI=True\n")

        outfile.write("\n")
        outfile.write(self.get_script_template(process_template))
        outfile.close()


def main(args = sys.argv):
    """This is the main body of the program that by default uses the arguments
from the command line."""

    parser = parser_setup()
    poptions = parser.parse_args()

    if poptions.quiet:
        logging.basicConfig(level=logging.WARNING, format=log_format)
    elif poptions.debug:
        logging.basicConfig(level=logging.DEBUG, format=log_format)
    else:
        # Set up the logging levels
        logging.basicConfig(level=logging.INFO, format=log_format)
        logging.getLogger("requests").setLevel(logging.WARNING)

    if not poptions.base_api_url and "LIMS_API_URL" in os.environ:
        api_url = os.environ["LIMS_API_URL"]
    elif poptions.base_api_url:
        api_url = poptions.base_api_url
    else:
        logging.error("Could not find LIMS API URL.\n")
        sys.exit(1)

    if not poptions.token and "LIMS_API_TOKEN" in os.environ:
        token = os.environ["LIMS_API_TOKEN"]
    elif poptions.token:
        token = poptions.token
    else:
        logging.error("Could not find LIMS API TOKEN.\n")
        sys.exit(1)

    process = ProcessSetUp(poptions, api_url, token)

    for align_id in poptions.align_ids:
        process.setup_alignment(align_id)

    if poptions.flowcell_label:
        process.setup_flowcell(poptions.flowcell_label)

    if poptions.tag:
        process.setup_tag(poptions.tag)

# This is the main body of the program that only runs when running this script
# doesn't run when imported, so you can use the functions above in the shell after importing
# without automatically running it
if __name__ == "__main__":
    main()