import base64
import re
import gzip
import file_io
import impute
import prs
import sys
import os
import shutil
import logging
from datetime import datetime, timezone

logger = logging.getLogger("app_logger")

logger.setLevel(logging.DEBUG)  # or INFO, WARNING, etc.
if not logger.handlers:
    handler = logging.StreamHandler()
    formatter = logging.Formatter(
        '[%(asctime)s] [%(levelname)s] %(message)s'
    )
    handler.setFormatter(formatter)
    logger.addHandler(handler)

#Extract the body from an event and strip JSON/File meta data
def extract(event):
    if 'body' in event:
        body = event['body']

        #Check if base64 encoded and try to decompress. If not it fails and we just decode
        if event['isBase64Encoded']:
                body = base64.b64decode(body)
        try:
            body = gzip.decompress(body).decode('utf-8')
            #return "decompress check"
        except gzip.BadGzipFile:
            pass
            body = body.decode('utf-8')

        body_lines = re.split(r'[\r\n]+', body)
        body_lines_filter = [element for element in body_lines if element.startswith('rs') or element.startswith('i')]
    else:
        body_lines_filter = []


    return body_lines_filter

#Convert each line in a list from csv to tsv. 
#Only use this function on lines already pre-processed with extract()
def convert_to_tsv(body):
    converted = []
    for line in body:
        # Ensure it's a data line and not malformed
        if ',' in line:
            converted.append(line.replace(',', '\t'))
        else:
            converted.append(line)  # already tab-separated (or something else)

    return converted

def guess_file_format(body):
    if not body:
        logger.error(f"[:Error:] Uploaded file is empty.")
        return 'error'
    line = body[0]
    columns = line.split('\t')
    nc = len(columns)
    if  nc == 4:
        return '23andme'
    elif nc == 5:
        return 'ancestry'
    else:
        logger.error(f"[:Error:] Unexpected number of columns in file: {nc}")
        return 'error'

#Count the number of variants in a file. 
def count_variants(vcf_file_path):
    variant_count = 0
    is_gzipped = vcf_file_path.endswith('.gz')
    open_func = gzip.open if is_gzipped else open
    with open_func(vcf_file_path, 'rt', encoding='utf-8') as file:
        for line in file:
            if not line.startswith('#'):
                variant_count += 1
    return variant_count

def clean_up(dir):
    if not os.path.exists(dir):
        logger.error(f"Directory '{dir}' does not exist.")
        return
    
    for item in os.listdir(dir):
        item_path = os.path.join(dir, item)

        try:
            if os.path.isfile(item_path):
                os.remove(item_path)  # Remove the file
                logger.debug(f"[DEBUG]: Removed file: {item_path}")
        except Exception as e:
            logger.error(f"Error removing {item_path}: {e}")


def handler(event, context):
    #clean_up('/tmp/')
    faipath = "/mnt/ref/ref/human_g1k_v37.fasta.fai"
    fapath = "/mnt/ref/ref/human_g1k_v37.fasta"
    vcfRef = "/mnt/ref/ref/1kgreference.bcf"
    mapFile = "/mnt/ref/ref/genetic_map_hg19_withX.txt.gz"
    haplo_ref_suffix = '1000g.Phase3.v5.With.Parameter.Estimates.msav'
    #Step 1: extract the uploaded file from the event
    body = extract(event)
    #Step2: Convert csv to tsv or keep tsv
    body = convert_to_tsv(body)
    #Step3: Convert to VCF

##TODO: [ERROR] IndexError: list index out of range
#Traceback (most recent call last):
#  File "/var/task/lambda.py", line 111, in handler
#    filetype = guess_file_format(body)
#  File "/var/task/lambda.py", line 61, in guess_file_format
#    line = body[0]

    filetype = guess_file_format(body)
    logger.debug(f"[DEBUG]: Version 0.2a")
    logger.debug(f"[DEBUG]: Guessed file format: {filetype}")
    row_count = len(body)
    logger.debug(f"[DEBUG]: Called Handler. Received object: {row_count} rows")
    if(filetype == '23andme'):
        snps = file_io.load_23andme_data(body)
    elif(filetype == 'ancestry'):
        snps = file_io.load_ancestry_data(body)
    else:
         logger.error(f"[ERROR] Failed to load genotypes. Exit.")
         sys.exit(1)

    fai = file_io.load_fai(faipath)
    records = file_io.get_vcf_records(snps, fai, fapath)
    infile = '/tmp/input.vcf'
    file_io.write_vcf(infile, records)
    row_count_vcf = file_io.count_vcf(infile)
    logger.debug(f"[DEBUG]: File conversion complete. VCF has {row_count_vcf} rows")
    chr = str(impute.extract_chromosome(infile))
    logger.debug(f"[DEBUG]: Chromosome found: {chr}")
    impute.index_vcf(infile)

    infile = '/tmp/input.vcf.gz'
    fileroot = '/mnt/ref/ref/'
    #file_io.list_files_recursive(fileroot)
    # Step 4: Pre-phasing
    logger.debug(f"[DEBUG]: Start phasing")
    impute.prePhase(infile, vcfRef, mapFile, chr)

    infile = '/tmp/phased.vcf.gz'
    impute.index_vcf(infile)

    #Step 3: Impute
    logger.debug(f"[DEBUG]: Start imputing")
    impute.impute(infile, haplo_ref_suffix, chr)

    #Step 4: Calculate Score
    infile = '/tmp/imputed.vcf.gz'
    row_count_imputed_vcf = file_io.count_vcf(infile)
    logger.debug(f"[DEBUG]: Imputing complete. Imputed VCF has {row_count_imputed_vcf} variants")
    prs_chr = prs.calc(infile, fileroot, chr)
    logger.debug(f"[DEBUG]: Calculated PRS for chr{chr}: {prs_chr}")
    #date_prefix = datetime.now(timezone.utc).strftime("%Y-%m-%d")
    #url = file_io.upload_file_to_s3(
    #    bucket_name="prs-tool",
    #    s3_key=f"prs_tool_debug/vcf/{date_prefix}_chr{chr}.vcf",
    #    local_file_path="/tmp/imputed.vcf.gz"
    #)
    #logger.debug(f"[DEBUG]: Stored imputed vcf as {url}")
    logger.debug(f"[DEBUG]: Cleaning up...")
    clean_up('/tmp/')
    logger.debug(f"[DEBUG]:All complete. Returning {list(prs_chr.keys())}")


    return file_io.dump(prs_chr, indent=2)
