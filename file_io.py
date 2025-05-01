import base64
import json
import boto3
import os
import numpy as np
import logging

#from pathlib import Path

logger = logging.getLogger("app_logger")


from botocore.exceptions import ClientError
# File I/O package to read&write various data objects to/from disk.

def list_files_recursive(directory):
    for root, dirs, files in os.walk(directory):
        for file in files:
            print(os.path.join(root, file))

def load_fai(faipath):

    index = {}
    with open(faipath) as f:
        for line in f:
            toks = line.split('\t')
            chrom = toks[0]
            if chrom == 'MT':
                chrom = 'M'
            length = int(toks[1])
            start = int(toks[2])
            linebases = int(toks[3])
            linewidth = int(toks[4])
            index[chrom] = (start, length, linebases, linewidth)
        return index

#Write body to file
def write_file(body, file_path):
    with open(file_path, 'w') as file:
        file.write(body)
    return file_path

#Read a file
def read_file(file_path, gzipped):
    if gzipped:
        with gzip.open(file_path, 'rt') as file:
            content = file.read()
    else:
        with file_path as file:
            content = file.read()
    return content

#Function to enable handling of ndarrays 
def ndarray_to_list(obj):
    if isinstance(obj, np.ndarray):
        return obj.tolist()
    raise TypeError(f"Object of type {type(obj)} is not JSON serializable")

#Dumps content into a json for testing
def dump(content, indent=2):

    return {
        'statusCode': 200,
        'headers': {
            "Access-Control-Allow-Origin": "*",
            "Content-Type": "application/json",
            "Access-Control-Allow-Methods": "GET,PUT,POST,DELETE,PATCH,OPTIONS",
            "Access-Control-Allow-Headers": "Content-Type,X-Amz-Date,Authorization,X-Api-Key,X-Amz-Security-Token"
    },
        'body': json.dumps(content, default=ndarray_to_list, indent=indent)
    }

def write_vcf_header(f):
    f.write(
"""##fileformat=VCFv4.2
##source=23andme_ancestryDNA_to_vcf
##reference=GRCh37
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
""")

#Write VCF to disk, so it can be read by
def write_vcf(outfile, records):
    snps = set()
    with open(outfile, 'w') as f:
        write_vcf_header(f)
        for record in records:
            if not (record[0] + ':' + record[1]) in snps:
                snps.add(record[0] + ':' + record[1]) # Skip duplicated
                f.write('\t'.join(record) + '\n')

## We take the loaded snplist from RAM (Ancestry)
def load_ancestry_data(lines):
    for line in lines:
        if line.startswith('#'): continue
        if line.strip():
            rsid, chrom, pos, allele1, allele2 = line.strip().split('\t')
            genotype = allele1 + allele2
            if chrom == 'MT':
                chrom = 'M'
            if genotype != '--':
                skip = False
                for x in genotype:
                    if x not in 'ACTG':
                        skip = True
                if not skip:
                    yield rsid, chrom, int(pos) - 1, genotype # subtract one because positions are 1-based indices

## We take the loaded snplist from RAM (23andMe)
def load_23andme_data(lines):
        for line in lines:
            if line.startswith('#'): continue
            if line.strip():
                rsid, chrom, pos, genotype = line.strip().split('\t')
                if chrom == 'MT':
                    chrom = 'M'
                if genotype != '--':
                    skip = False
                    for x in genotype:
                        if x not in 'ACTG':
                            skip = True
                    if not skip:
                        yield rsid, chrom, int(pos) - 1, genotype # subtract one because positions are 1-based indices

def check_efs_file(file_path):

    if os.path.exists(file_path):
        rows = []
        with open(file_path, 'r') as file:
            for i in range(10):
                line = file.readline()
                if not line:
                    break
                rows.append(line.strip())
            return rows
    else:
        logger.error(f"[Error] {file_path} not found!")
        raise FileNotFoundError()
        #return 'Error: File not found!'

def get_vcf_records(pos_list, fai, fapath):

    #Function returns the alt allele aligned to +
    def get_alts(ref, genotype):
            for x in genotype:
                assert x in 'ACGT'

            if len(genotype) == 1:
#                if ref in genotype:
#                    return [genotype]
                return [genotype]

            if ref == genotype[0] and ref == genotype[1]:
                return [genotype[0]] # we always geturn a genotype
            if ref == genotype[0]:
                return [genotype[1]]
            if ref == genotype[1]:
                return [genotype[0]]
            return [genotype[0], genotype[1]]

    # Iterate over each tuple in pos_list
    with open(fapath) as f:
        for (rsid, chrom, pos, genotype) in pos_list:
            start, _, linebases, linewidth = fai[chrom]
            n_lines = int(pos / linebases)
            n_bases = pos % linebases
            n_bytes = start + n_lines * linewidth + n_bases
            f.seek(n_bytes)
    # Stream the reference allele from the FASTA file
            ref = f.read(1)
            alts = get_alts(ref, genotype)
            pos = str(pos + 1)
            diploid = len(genotype) == 2
#            assert ref not in alts # we always geturn a genotype
            assert len(alts) <= 2
            if diploid:
                if len(alts) == 2:
                #Drop multi alleles
                    if alts[0] == alts[1]:
                        yield (chrom, pos, rsid, ref, alts[0], '.', '.', '.', 'GT', '1/1')
                elif len(alts) == 1:
                    yield (chrom, pos, rsid, ref, alts[0], '.', '.', '.', 'GT', '0/1')
            elif len(alts) == 1:
                yield (chrom, pos, rsid, ref, alts[0], '.', '.', '.', 'GT', '1')

def count_vcf(file_path: str):
    import gzip
    open_fn = gzip.open if file_path.endswith('.gz') else open
    variant_count = 0
    with open_fn(file_path, 'rt') as file:
        for line in file:
            # Skip comments and headers
            if line.startswith('#'):
                continue
            # Count the variant lines
            variant_count += 1

    return variant_count

def upload_file_to_s3(bucket_name: str, s3_key: str, local_file_path: str) -> bool:
     
    s3 = boto3.client('s3')
     # Upload the file
    logger.debug(f"[:DEBUG:]Uploading {local_file_path} to s3://{bucket_name}/{s3_key}")
    s3.upload_file(
        Filename=local_file_path,
        Bucket=bucket_name,
        Key=s3_key,
        ExtraArgs={'ContentType': 'text/plain'}  # Adjust if needed
    )
    logger.debug(f"[:DEBUG:]Uploading successful!")
    # Generate a presigned URL (valid for 1 hour)
    presigned_url = s3.generate_presigned_url(
        ClientMethod='get_object',
        Params={'Bucket': bucket_name, 'Key': s3_key},
        ExpiresIn=3600
    )
    return presigned_url
