import subprocess
import os
import logging
import numpy as np

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
        logger.error(f"No chromosome found in the VCF file.")
        raise ChromosomeCountError()

    return chromosomes.pop()

def index_vcf(vcf_path):

    #File check
    if not os.path.isfile(vcf_path):
        logger.error(f"The file {vcf_path} does not exist.")
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
        logger.error(f"An error occurred while running the command: {e.cmd}")
        logger.error(f"{e.stderr}")
        raise
    except Exception as e:
        logger.error(f"An unexpected error occurred: {str(e)}")
        raise
