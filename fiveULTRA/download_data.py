import requests
import zipfile
import argparse
import logging
import hashlib
import shutil
import time
import sys
from pathlib import Path

# --- Configuration ---
# URL for the data download
DATA_URL = 'https://hgidsoft.rockefeller.edu/5ULTRA/5ultra_data_archive.zip'
# Number of times to retry a failed download
DOWNLOAD_RETRIES = 3
# Seconds to wait between retries
RETRY_DELAY = 5

def setup_logging():
    """Configures logging to print messages to the console."""
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        stream=sys.stdout,
    )

def calculate_sha256(filepath: Path) -> str:
    """Calculates and returns the SHA256 hash of a file."""
    sha256_hash = hashlib.sha256()
    with open(filepath, "rb") as f:
        # Read and update hash in chunks of 4K
        for byte_block in iter(lambda: f.read(4096), b""):
            sha256_hash.update(byte_block)
    return sha256_hash.hexdigest()

def main():
    """
    Downloads, verifies, and extracts the necessary data for the 5ULTRA tool.
    """
    parser = argparse.ArgumentParser(
        description="Download and set up necessary data for 5ULTRA.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    default_data_dir = Path.home() / ".5ULTRA" / "data"
    parser.add_argument(
        '--data-dir',
        type=Path,
        default=default_data_dir,
        help='Path to the directory where data will be stored.'
    )
    parser.add_argument(
        '--force',
        action='store_true',
        help='Force re-download and extraction even if the data directory already exists.'
    )
    parser.add_argument(
        '--no-cleanup',
        action='store_true',
        help='Do not delete the downloaded ZIP file after extraction.'
    )
    args = parser.parse_args()

    data_dir: Path = args.data_dir
    local_zip_filename = 'data.zip'
    local_zip_filepath = data_dir / local_zip_filename

    # Check if data directory already exists and is not empty
    if data_dir.exists() and any(data_dir.iterdir()) and not args.force:
        logging.info(f"Data directory '{data_dir}' already exists and is not empty.")
        logging.info("Use --force to re-download and overwrite.")
        sys.exit(0)

    # Ensure the target directory exists
    data_dir.mkdir(parents=True, exist_ok=True)

    # --- Download Phase ---
    logging.info(f"Attempting to download data from {DATA_URL}")
    for attempt in range(DOWNLOAD_RETRIES):
        try:
            with requests.get(DATA_URL, stream=True, timeout=30) as r:
                r.raise_for_status()
                total_size = int(r.headers.get('content-length', 0))
                
                # Import tqdm locally to keep dependency optional for non-download scenarios
                from tqdm import tqdm
                
                logging.info(f"Downloading to {local_zip_filepath}...")
                with open(local_zip_filepath, 'wb') as f, tqdm(
                    desc="Downloading",
                    total=total_size,
                    unit='iB',
                    unit_scale=True,
                    unit_divisor=1024,
                ) as bar:
                    for chunk in r.iter_content(chunk_size=8192):
                        size = f.write(chunk)
                        bar.update(size)
            logging.info("Download complete.")
            break  # Exit retry loop on success
        except requests.exceptions.RequestException as e:
            logging.error(f"Download attempt {attempt + 1} failed: {e}")
            if attempt + 1 < DOWNLOAD_RETRIES:
                logging.info(f"Retrying in {RETRY_DELAY} seconds...")
                time.sleep(RETRY_DELAY)
            else:
                logging.error("All download attempts failed. Please check your network connection and the URL.")
                sys.exit(1)

    # --- Verification Phase ---
    logging.info("Verifying file integrity...")
    try:
        zip_hash = calculate_sha256(local_zip_filepath)
        logging.info(f"SHA256 Hash of downloaded file: {zip_hash}")
    except IOError as e:
        logging.error(f"Could not read downloaded file for verification: {e}")
        sys.exit(1)


    # --- Extraction Phase (Transactional) ---
    temp_extract_dir = data_dir.with_name(data_dir.name + "_tmp")
    
    try:
        logging.info(f"Extracting data to a temporary location...")
        if temp_extract_dir.exists():
            shutil.rmtree(temp_extract_dir) # Clean up previous failed attempts
        temp_extract_dir.mkdir()

        with zipfile.ZipFile(local_zip_filepath, 'r') as zip_ref:
            zip_ref.extractall(temp_extract_dir)
        
        logging.info("Extraction successful. Moving data to final location.")
        # Remove old data directory if it exists
        if data_dir.exists():
            shutil.rmtree(data_dir)
        # Move new data into place
        shutil.move(str(temp_extract_dir), str(data_dir))

    except zipfile.BadZipFile:
        logging.error("Error: The downloaded file is not a valid ZIP archive.")
        shutil.rmtree(temp_extract_dir, ignore_errors=True) # Clean up temp dir
        sys.exit(1)
    except Exception as e:
        logging.error(f"An unexpected error occurred during extraction: {e}")
        shutil.rmtree(temp_extract_dir, ignore_errors=True) # Clean up temp dir
        sys.exit(1)
    
    # --- Cleanup Phase ---
    if not args.no_cleanup:
        logging.info("Cleaning up downloaded ZIP file...")
        try:
            local_zip_filepath.unlink()
        except OSError as e:
            logging.warning(f"Could not remove ZIP file: {e}")
    else:
        logging.info("Skipping cleanup of ZIP file as requested.")

    logging.info(f"Setup complete. Data is ready in '{data_dir}'.")


if __name__ == '__main__':
    setup_logging()
    try:
        main()
    except Exception as e:
        logging.critical(f"An unhandled error occurred: {e}", exc_info=True)
        sys.exit(1)