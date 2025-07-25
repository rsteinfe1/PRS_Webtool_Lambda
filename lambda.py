
import base64
import re
import gzip
import file_io
import impute
import prs
import sys
import os
import logging
import json
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

#Extract the genotyping from an event and strip JSON/File meta data
def extract(event):
    body_raw = event['body']
    if event.get('isBase64Encoded', False):
        body_raw = base64.b64decode(body_raw)

    body = json.loads(body_raw)
    return body

def getGTs(body:dict) -> list:
    """
    Extracts the 'genotypes' field from the body of the event.
    If 'genotypes' is not present, it returns an empty list.
    """
    if 'genotypes' in body:
        genotypes = body.get('genotypes', '')
        gt_lines = re.split(r'[\r\n]+', genotypes)
        gtlines_filter = [line for line in gt_lines if line.startswith('rs') or line.startswith('i')]
        return gtlines_filter
    else:
        logger.error(f"[:Error:] 'genotypes' field not found in the body.")
        return []

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
    clean_up('/tmp/')
    #logger.debug(f"[DEBUG]: Received event: {json.dumps(event)}")

    method = event.get("httpMethod")  # REST API
    if not method and "requestContext" in event:
        method = event["requestContext"].get("http", {}).get("method")  # HTTP API

    if method == "OPTIONS":
        return {
            "statusCode": 200,
            "headers": {
                "Access-Control-Allow-Origin": "*",
                "Access-Control-Allow-Methods": "POST,OPTIONS",
                "Access-Control-Allow-Headers": "Content-Type"
            },
            "body": ""
        }
    
    fai36path = "/mnt/ref/ref/human_genome_v36.fa.fai"
    fa36path = "/mnt/ref/ref/human_genome_v36.fa"
    chain1 = "/mnt/ref/ref/hg18ToHg19.over.chain.gz"
    chain2 = "/mnt/ref/ref/hg19ToHg38.over.chain.gz"
    faipath = "/mnt/ref/ref/human_g1k_v37.fasta.fai"
    fapath = "/mnt/ref/ref/human_g1k_v37.fasta"
    vcfRef = "/mnt/ref/ref/1kgreference.bcf"
    mapFile = "/mnt/ref/ref/genetic_map_hg19_withX.txt.gz"
    haplo_ref_suffix = '1000g.Phase3.v5.With.Parameter.Estimates.msav'
    logger.debug(f"[DEBUG]: Version 0.4a")
    #logger.debug(f"[DEBUG]: Received event: {json.dumps(event)}")
    #Step 1: extract the uploaded payload from the event
    body = extract(event)
    build = body.get('build', 'NA')  # Default to NA if not specified
    logger.debug(f"[DEBUG]: Received build: {build}")
    body = getGTs(body)
    logger.debug(f"[DEBUG]: Extracted genotypes: {len(body)} lines. Class {type(body)}. First lines: {body[:5]}")
    #Step2: Convert csv to tsv or keep tsv
    body = convert_to_tsv(body)
    #Step3: Convert to VCF

    filetype = guess_file_format(body)
    logger.debug(f"[DEBUG]: Guessed file format: {filetype}")
    row_count = len(body)
    logger.debug(f"[DEBUG]: Called file handler. Received object: {row_count} rows")
    if filetype == '23andme':
        snps = file_io.load_23andme_data(body)
    elif filetype == 'ancestry':
        snps = file_io.load_ancestry_data(body)
    else:
        #Todo: Write a build guesser function here.
        logger.error(f"[ERROR] Failed to load genotypes. Exit.")
        return {
                    'statusCode': 400,
                    'body': json.dumps({'Error': 'Unknown file format.'})
                }

    infile = '/tmp/input.vcf'
    if build == 'GRCh36':
        fai = file_io.load_fai(fai36path)
        records = file_io.get_vcf_records(snps, fai, fa36path)
        file_io.write_vcf(infile, records)
        impute.liftOver(chain1, infile, fapath)
    elif build == 'GRCh37':
        fai = file_io.load_fai(faipath)
        records = file_io.get_vcf_records(snps, fai, fapath)
        file_io.write_vcf(infile, records)
    else:
        #We determine build by matching rsID to Locus, lift over (if necessary) and write.
        locusid = f"/mnt/ref/ref/dbSNP_151_idlocus_hg19_chr{chr}.txt"
        matchbuild37 = impute.match_locusids_from_body(body, locusid)
        if matchbuild37 > 0.9:
            build = 'GRCh37'
            logger.debug(f"[DEBUG]: Detected build GRCh37 with {matchbuild37} match ratio.")
            fai = file_io.load_fai(faipath)
            records = file_io.get_vcf_records(snps, fai, fapath)
            file_io.write_vcf(infile, records)
        else:
            locusid = f"/mnt/ref/ref/dbSNP_151_idlocus_hg18_chr{chr}.txt"
            matchbuild36 = impute.match_locusids_from_body(body, locusid)
            if matchbuild36 > 0.9:
                logger.debug(f"[DEBUG]: Detected build GRCh37 with {matchbuild36} match ratio.")
                fai = file_io.load_fai(faipath)
                records = file_io.get_vcf_records(snps, fai, fapath)
                file_io.write_vcf(infile, records)
                impute.liftOver(chain1, infile, fapath)
            else:
                locusid = f"/mnt/ref/ref/dbSNP_151_idlocus_hg38_chr{chr}.txt"
                matchbuild38 = impute.match_locusids_from_body(body, locusid)
                if matchbuild38 > 0.9:
                    logger.debug(f"[DEBUG]: Detected build GRCh38 with {matchbuild38} match ratio.")
                    fai = file_io.load_fai(faipath)
                    records = file_io.get_vcf_records(snps, fai, fapath)
                    file_io.write_vcf(infile, records)
                    impute.liftOver(chain2, infile, fapath)
                else:
                    #Unknown build
                    logger.error(f"[ERROR] Failed to detect build from genotypes. Match ratios: GRCh36: {matchbuild36}, GRCh37: {matchbuild37}, GRCh38: {matchbuild38} ")
                    logger.error(f"[ERROR] Unknown build: {build}. Exiting.")
                    return {
                        'statusCode': 400,
                        'body': json.dumps({'Error': 'Unknown build'})
                    }

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
