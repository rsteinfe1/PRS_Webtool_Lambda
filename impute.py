import subprocess
import os
import logging
import sys
import numpy as np
import shutil
import gzip

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
    """Exception for incorrect number of chromosomes in file."""
    pass

class ChromosomeValueError(Exception):
    """Exception for wrong value of chromosomes in file."""
    pass

def extract_chromosome_from_vcf(vcf_file_path):
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

def extract_chromosome_from_body(body):  
    chromosomes = set()
    for line in body:
        if line.startswith("#") or not line.strip():
            continue
        fields = line.strip().split("\t")
        if len(fields) >= 2:
            chrom = fields[1].strip()
            if chrom:
                chromosomes.add(chrom)

    if len(chromosomes) != 1:
        logger.error(f"Expected exactly one chromosome, found {len(chromosomes)}: {chromosomes}")
        raise ChromosomeCountError()

    chrom = chromosomes.pop()
    if not chrom.isdigit() or not (1 <= int(chrom) <= 22):
        logger.error(f"Chromosome '{chrom}' is not a valid chromosome (1-22).")
        raise ChromosomeValueError()
    
    return chrom


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
        
        return bgzipped_file
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
    Adds missing ##contig and ##FORMAT=GT lines to the VCF header.
    Avoids duplicate entries. Modifies the VCF file in place.
    """
    open_func = gzip.open if vcf_file.endswith(".gz") else open
    write_mode = "wt" if vcf_file.endswith(".gz") else "w"

    # Read contigs from .fai file
    contig_lines = []
    with open(fai_file, "r") as fai:
        for line in fai:
            chrom, length = line.strip().split("\t")[:2]
            contig_lines.append(f"##contig=<ID={chrom},length={length}>\n")

    contig_set = set(cl.strip() for cl in contig_lines)
    gt_format_line = '##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">\n'

    header = []
    body = []
    existing_contigs = set()
    has_gt = False

    with open_func(vcf_file, "rt") as vcf:
        for line in vcf:
            if line.startswith("##contig="):
                existing_contigs.add(line.strip())
            elif line.startswith("##FORMAT=<ID=GT"):
                has_gt = True

            if line.startswith("#CHROM"):
                # Inject missing contig and GT FORMAT lines before #CHROM
                missing_contigs = [cl for cl in contig_lines if cl.strip() not in existing_contigs]
                header.extend(missing_contigs)
                if not has_gt:
                    header.append(gt_format_line)
                header.append(line)
                break
            else:
                header.append(line)

        body = list(vcf)  # Read remaining lines (body)

    logger.debug(f"Injected contigs (added {len(missing_contigs)}) to {vcf_file}")
    if not has_gt:
        logger.debug(f"Injected GT FORMAT line to {vcf_file}")

    # Write updated VCF
    with open_func(vcf_file, write_mode) as out:
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
    run_cmd(f"/usr/local/bcftools-1.22/bcftools sort {lifted_raw} -Ov -o {sorted_tmp}", shell=True)

    # Step 4: Rename to final path in /tmp/ before compression
    shutil.move(sorted_tmp, input_vcf)

def match_locusids_from_body(body, locusid_file):

    # Step 1: Validate and parse body
    locus_keys = set()
    for line in body:
        if line.startswith("#") or not line.strip():
            continue
        fields = line.strip().split("\t")
        if len(fields) < 3:
            continue  # skip malformed lines
        rsid, chrom, pos = fields[0], fields[1], fields[2]
        if not rsid or not chrom or not pos or not rsid.startswith('rs'):
            continue
        locus_keys.add(f"{rsid}:{chrom}:{pos}")

    total = len(locus_keys)
    if total == 0:
        logger.warning("Could not parse body, body is empty or no dbSNP identifiers found.")
        return 0.0

    # Step 2: Check if file exists and is readable
    if not os.path.exists(locusid_file):
        logger.error(f"LocusID file {locusid_file} not found.")
        return 0.0
    if not os.path.isfile(locusid_file):
        logger.error(f"{locusid_file} is not a regular file.")
        return 0.0

    try:
        with open(locusid_file, "r") as f:
            valid_locusids = {line.strip() for line in f if line.strip()}
    except Exception as e:
        logger.error(f" Failed reading {locusid_file}: {e}")
        return 0.0

    # Step 3: Compute match ratio
    matches = locus_keys & valid_locusids
    ratio = len(matches) / total if total > 0 else 0.0

    return ratio

def print_vcf_preview(vcf_path, n=10, show_header=True):
    """
    Prints the first `n` data rows from a VCF file.
    
    Parameters:
        vcf_path (str): Path to .vcf or .vcf.gz file.
        n (int): Number of data lines (non-header) to print.
        show_header (bool): If True, also print header lines.
    """
    open_func = gzip.open if vcf_path.endswith(".gz") else open

    count = 0
    with open_func(vcf_path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                if show_header:
                    print(line.strip())
                continue
            print(line.strip())
            count += 1
            if count >= n:
                break

    if count == 0:
        print("No data lines found in the VCF.")


def normalize_vcf(vcf_file, fa_file, fai_file):
    logger.debug(f"Normalizing: {vcf_file}")
    
    try:
        # Inject contigs from the reference FASTA index
        inject_contigs(vcf_file, fai_file)
        # Match alleles and sort the VCF file using bcftools
        norm_vcf = '/tmp/normalized.sorted.vcf.gz'
        run_cmd(f"/usr/local/bcftools-1.22/bcftools -m -both --threads 4 -f {fa_file} -cs -Ov {vcf_file} | \
                /usr/local/bcftools-1.22/bcftools +fixref {vcf_file} --threads 4 -Ov -- -f {fa_file} -m swap | \
                /usr/local/bcftools-1.22/bcftools sort -Oz -o {norm_vcf}", shell=True)

        # Index the sorted VCF
        index_vcf(norm_vcf)
        return norm_vcf
    
    except Exception as e:
        logger.error(f"Normalizing VCF failed: {str(e)}")
        raise

