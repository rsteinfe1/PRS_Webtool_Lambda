import subprocess
import os
import logging
import sys
import numpy as np
import shutil

logger = logging.getLogger("app_logger")

#Function to run phasing with Eagle
def prePhase(vcfInput, vcfRef, mapFile, chrom):

    command = ['eagle', '--vcfRef', vcfRef,
               '--vcfTarget', vcfInput,
               '--geneticMapFile', mapFile,
               '--outPrefix', '/tmp/phased',
               '--allowRefAltSwap',
                '--vcfOutFormat', 'z',
                '--numThreads', '10',
                '--chrom', chrom ]
    try:
        subprocess.run(command, text=True, check=True)
    except subprocess.CalledProcessError as e:
        logger.error(f"[ERROR] Command '{e.cmd}' returned non-zero exit status {e.returncode}")
        logger.error(f"Error output: {e.stderr}")
        raise
    except Exception as e:
        logger.error(f"An exception occurred: {str(e)}")

def impute(vcfInput, haplo_ref_suffix, chr):
    command = [ 'minimac4', '--output', "/tmp/imputed.vcf.gz",
                '--threads', '10',
                '--format', 'GT,DS,GP',
                '--all-typed-sites',
                '--empirical-output', f"/tmp/empiricalDosage.vcf.gz",
                f"/mnt/ref/ref/{chr}.{haplo_ref_suffix}", vcfInput ]
    try:
        subprocess.run(command, text=True, check=True)

    except subprocess.CalledProcessError as e:
         logger.error(f"Command '{e.cmd}' returned non-zero exit status {e.returncode}")
         logger.error(f"Error output: {e.stderr}")
    except Exception as e:
         logger.error(f"An exception occurred: {str(e)}")

class ChromosomeCountError(Exception):
    """Exception for incorrect number of chromosomes in VCF file."""
    pass

def extract_chromosome(vcf_file_path):
    chromosomes = set()
    with open(vcf_file_path, 'r') as file:
        for line in file:
            if line.startswith('#'):
                continue
            columns = line.split('\t')
            chromosome = columns[0]
            chromosomes.add(chromosome)
            chrom_len = len(chromosomes)
            if chrom_len > 1:
                logger.error(f"More than one chromosome ({chrom_len}) found in the VCF file: {chromosomes}")
                raise ChromosomeCountError()

    if len(chromosomes) == 0:
        logger.error(f"[ERROR] No chromosome found in the VCF file.")
        raise ChromosomeCountError()

    return chromosomes.pop()

def index_vcf(vcf_path):

    #File check
    if not os.path.isfile(vcf_path):
        logger.error(f"[ERROR] The file {vcf_path} does not exist.")
        raise FileNotFoundError()

    is_bgzipped = vcf_path.endswith('.gz')
    bgzipped_file = vcf_path if is_bgzipped else vcf_path + '.gz'

    try:
        if not is_bgzipped:
            bgzip_command = ['bgzip', '-c', vcf_path]
            with open(bgzipped_file, 'wb') as f_out:
                subprocess.run(bgzip_command, stdout=f_out, check=True)
        # Check if tabix index already exists
        tabix_file = bgzipped_file + '.tbi'
        if not os.path.isfile(tabix_file):
            tabix_command = ['tabix', '-fp', 'vcf', bgzipped_file]
            subprocess.run(tabix_command, check=True)
    except subprocess.CalledProcessError as e:
        logger.error(f"[ERROR] Failed running: {e.cmd}")
        logger.error(f"{e.stderr}")
        raise
    except Exception as e:
        logger.error(f"[ERROR] An unexpected error occurred: {str(e)}")
        raise

def run_cmd(cmd, shell=False):
    print(f"Running: {cmd if isinstance(cmd, str) else ' '.join(cmd)}")
    result = subprocess.run(cmd, shell=shell, capture_output=True, text=True)
    if result.returncode != 0:
        logger.error(f"[ERROR] Failed running: {result.stderr}")
        sys.exit(1)

def inject_contigs(vcf_file, fai_file):
    """
    Injects contig headers from a .fai file into a VCF file.
    This is necessary to ensure that the VCF file can be sorted and indexed correctly.
    print(f"Running: {cmd if isinstance(cmd, str) else ' '.join(cmd)}")
    """
    contig_lines = []
    with open(fai_file) as fai:
        for line in fai:
            chrom, length = line.strip().split("\t")[:2]
            contig_lines.append(f"##contig=<ID={chrom},length={length}>\n")

    header = []
    body = []
    with open(vcf_file) as vcf:
        for line in vcf:
            if line.startswith("##"):
                header.append(line)
            elif line.startswith("#CHROM"):
                header.extend(contig_lines)
                header.append(line)
                break
            else:
                header.append(line)
        body = vcf.readlines()

    with open(vcf_file, "w") as out:
        out.writelines(header)
        out.writelines(body)


def liftOver(chain, input_vcf, ref_fasta):
    # File prefix for safe naming. 
    prefix = os.path.splitext(os.path.basename(input_vcf))[0].replace('.vcf', '')

    # Temporary working files in /tmp
    lifted_raw = f"/tmp/{prefix}.lifted.unsorted.vcf"
    sorted_tmp = f"/tmp/{prefix}.sorted.tmp.vcf"
    
    # Step 1: CrossMap
    run_cmd(["CrossMap", "vcf", chain, input_vcf, ref_fasta, lifted_raw])

    # Step 2: Inject correct contig headers from .fai
    inject_contigs(lifted_raw, ref_fasta + ".fai")

    # Step 3: Sort with bcftools
    run_cmd(f"bcftools sort {lifted_raw} -o {sorted_tmp}", shell=True)

    # Step 4: Rename to final path in /tmp/ before compression
    shutil.move(sorted_tmp, input_vcf)




