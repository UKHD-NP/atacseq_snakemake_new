from urllib.request import urlretrieve
import gzip
import shutil
import ssl
import os
import sys


def info(msg):
    print(f"[INFO] {msg}")


def warning(msg):
    print(f"[WARNING] {msg}")


def fatal(msg):
    sys.exit(f"[ERROR] {msg}")


REFERENCE_URLS = {
    "hg19": {
        "fasta": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/GRCh37.p13.genome.fa.gz",
        "gtf": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_19/gencode.v19.annotation.gtf.gz",
        "blacklist": "https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/hg19-blacklist.v2.bed.gz",
    },
    "hg38": {
        "fasta": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/GRCh38.primary_assembly.genome.fa.gz",
        "gtf": "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_46/gencode.v46.primary_assembly.basic.annotation.gtf.gz",
        "blacklist": "https://raw.githubusercontent.com/Boyle-Lab/Blacklist/master/lists/hg38-blacklist.v2.bed.gz",
    },
}


assembly = config["ref"]["assembly"]
user_fasta = config["ref"].get("fasta", "")
user_gtf = config["ref"].get("gtf", "")
user_blacklist = config["ref"].get("blacklist", "")
url_fasta = None
url_gtf = None
url_blacklist = None


if assembly == "custom":
    if not user_fasta or not user_gtf or not user_blacklist:
        fatal(
            "For custom assembly configurations, please specify ref.fasta, ref.gtf, "
            "and ref.blacklist in config/config.yml."
        )
elif assembly in REFERENCE_URLS:
    url_fasta = REFERENCE_URLS[assembly]["fasta"]
    url_gtf = REFERENCE_URLS[assembly]["gtf"]
    url_blacklist = REFERENCE_URLS[assembly]["blacklist"]
else:
    fatal(f"Invalid 'ref.assembly': '{assembly}'. Allowed values: hg19, hg38, custom.")


# Keep references in a deterministic location inside the project.
ref_dir = os.path.join("references", assembly)
os.makedirs(ref_dir, exist_ok=True)


# Canonical reference paths used by downstream rules.
config["ref"]["fasta"] = os.path.join(ref_dir, config["ref"]["assembly"] + ".fa")
config["ref"]["gtf"] = os.path.join(ref_dir, config["ref"]["assembly"] + ".gtf")
config["ref"]["blacklist"] = os.path.join(ref_dir, config["ref"]["assembly"] + ".blacklist.bed")
config["ref"]["bed"] = os.path.join(ref_dir, config["ref"]["assembly"] + ".bed")
config["ref"]["chromsizes"] = os.path.join(ref_dir, config["ref"]["assembly"] + ".sizes")
config["ref"]["autosomes"] = os.path.join(ref_dir, config["ref"]["assembly"] + ".autosomes.txt")
config["ref"]["tss"] = os.path.join(ref_dir, config["ref"]["assembly"] + ".tss.bed")
config["ref"]["include_regions"] = os.path.join(ref_dir, config["ref"]["assembly"] + ".include_regions.bed")


# Helper function to decompress files
def decompress_file(input_file, output_file):
    """Decompress gzipped file to output path with error handling."""
    try:
        with gzip.open(input_file, "rb") as f_in:
            with open(output_file, "wb") as f_out:
                shutil.copyfileobj(f_in, f_out)
        return True
    except Exception as e:
        warning(f"Failed to decompress '{input_file}': {str(e)}")
        return False


# Helper function to copy or decompress files
def stage_local_reference(source_path, target_path, label):
    """Copy or decompress a local reference file into the references directory."""
    if not os.path.isfile(source_path):
        fatal(f"Custom {label} file does not exist: {source_path}")

    if os.path.abspath(source_path) == os.path.abspath(target_path):
        info(f"Using existing {label} file in place: {source_path}")
        return

    os.makedirs(os.path.dirname(target_path), exist_ok=True)
    
    if source_path.endswith(".gz"):
        info(f"Decompressing custom {label} file: {source_path}")
        if not decompress_file(source_path, target_path):
            fatal(f"Could not decompress custom {label} file: {source_path}")
    else:
        try:
            shutil.copyfile(source_path, target_path)
        except Exception as e:
            fatal(f"Could not copy custom {label} file '{source_path}' -> '{target_path}': {str(e)}")
            

# Helper function to download files
def download_reference_file(url, compressed_path, target_path, label):
    """Download and decompress a remote reference file."""
    try:
        info(f"Downloading {label} reference from: {url}")
        # Temporarily disable SSL verification only for this download call
        _orig_ctx = ssl._create_default_https_context
        ssl._create_default_https_context = ssl._create_unverified_context
        try:
            urlretrieve(url, compressed_path)
        finally:
            ssl._create_default_https_context = _orig_ctx
        info(f"Decompressing downloaded {label} file: {compressed_path}")
        if not decompress_file(compressed_path, target_path):
            fatal(f"Could not decompress downloaded {label} file: {compressed_path}")
        os.remove(compressed_path)
    except Exception as e:
        fatal(f"Could not download {label} file from '{url}': {str(e)}")


# Check if FASTA reference already exists. Download it if not.
if not os.path.isfile(config["ref"]["fasta"]):
    if config["ref"]["assembly"] == "custom":
        info("Using custom FASTA reference.")
        stage_local_reference(user_fasta, config["ref"]["fasta"], "FASTA")
    else:
        info(f"Downloading {config['ref']['assembly']} FASTA reference.")
        fasta_compressed = os.path.join(ref_dir, config["ref"]["assembly"] + ".fa.gz")
        download_reference_file(url_fasta, fasta_compressed, config["ref"]["fasta"], "FASTA")


# Check if GTF annotation already exists. Download it if not.
if not os.path.isfile(config["ref"]["gtf"]):
    if config["ref"]["assembly"] == "custom":
        info("Using custom GTF annotation.")
        stage_local_reference(user_gtf, config["ref"]["gtf"], "GTF")
    else:
        info(f"Downloading {config['ref']['assembly']} GTF annotation.")
        gtf_compressed = os.path.join(ref_dir, config["ref"]["assembly"] + ".gtf.gz")
        download_reference_file(url_gtf, gtf_compressed, config["ref"]["gtf"], "GTF")


# Check if blacklist already exists. Download it if not.
if not os.path.isfile(config["ref"]["blacklist"]):
    if config["ref"]["assembly"] == "custom":
        info("Using custom blacklist annotation.")
        stage_local_reference(user_blacklist, config["ref"]["blacklist"], "BLACKLIST")
    else:
        info(f"Downloading {config['ref']['assembly']} blacklist annotation.")
        blacklist_compressed = os.path.join(ref_dir, config["ref"]["assembly"] + ".blacklist.bed.gz")
        download_reference_file(url_blacklist, blacklist_compressed, config["ref"]["blacklist"], "BLACKLIST")
