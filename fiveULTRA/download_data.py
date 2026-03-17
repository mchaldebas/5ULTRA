import requests
import zipfile
import argparse
import logging
import hashlib
import shutil
import time
import os
import sys
from pathlib import Path

# --- Configuration ---
# Extracted from your link: https://drive.google.com/file/d/1zVYEZt72qgCgye8_Az6zfAmUk1-cGp9U/
GOOGLE_DRIVE_FILE_ID = '1zVYEZt72qgCgye8_Az6zfAmUk1-cGp9U'
DOWNLOAD_RETRIES = 3
RETRY_DELAY = 5

def setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        stream=sys.stdout,
    )

def get_confirm_token(response):
    """Helper to extract the Google Drive large file confirmation token."""
    for key, value in response.cookies.items():
        if key.startswith('download_warning'):
            return value
    return None

def calculate_sha256(filepath: Path) -> str:
    sha256_hash = hashlib.sha256()
    with open(filepath, "rb") as f:
        for byte_block in iter(lambda: f.read(4096), b""):
            sha256_hash.update(byte_block)
    return sha256_hash.hexdigest()

def download_from_google_drive(file_id, destination):
    """Downloads large files from Google Drive by handling the confirmation prompt."""
    URL = "https://docs.google.com/uc?export=download"
    session = requests.Session()
    
    # First request to check for the 'large file' warning
    response = session.get(URL, params={'id': file_id}, stream=True)
    token = get_confirm_token(response)

    if token:
        # Second request with the token to actually start the download
        params = {'id': file_id, 'confirm': token}
        response = session.get(URL, params=params, stream=True)
    
    return response

def main():
    parser = argparse.ArgumentParser(
        description="Download and set up necessary data for 5ULTRA from Google Drive.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    default_data_dir = Path.home() / ".5ULTRA" / "data"
    parser.add_argument('--data-dir', type=Path, default=default_data_dir)
    parser.add_argument('--force', action='store_true')
    parser.add_argument('--no-cleanup', action='store_true')
    args = parser.parse_args()

    data_dir: Path = args.data_dir
    local_zip_filename = 'data.zip'
    local_zip_filepath = data_dir / local_zip_filename

    if data_dir.exists() and any(data_dir.iterdir()) and not args.force:
        logging.info(f"Data directory '{data_dir}' already exists. Use --force to overwrite.")
        sys.exit(0)

    data_dir.mkdir(parents=True, exist_ok=True)

    # --- Download Phase ---
    logging.info(f"Connecting to Google Drive (ID: {GOOGLE_DRIVE_FILE_ID})...")
    
    success = False
    for attempt in range(DOWNLOAD_RETRIES):
        try:
            response = download_from_google_drive(GOOGLE_DRIVE_FILE_ID, local_zip_filepath)
            response.raise_for_status()
            
            # Use 'Content-Length' if available (G-Drive doesn't always send it for large files)
            total_size = int(response.headers.get('content-length', 0))
            
            from tqdm import tqdm
            
            logging.info(f"Downloading to {local_zip_filepath}...")
            with open(local_zip_filepath, 'wb') as f, tqdm(
                desc="Progress",
                total=total_size if total_size > 0 else None,
                unit='iB',
                unit_scale=True,
                unit_divisor=1024,
            ) as bar:
                for chunk in response.iter_content(chunk_size=32768): # Larger chunks for 3GB file
                    if chunk: 
                        f.write(chunk)
                        bar.update(len(chunk))
            
            logging.info("Download complete.")
            success = True
            break 
        except Exception as e:
            logging.error(f"Attempt {attempt + 1} failed: {e}")
            if attempt + 1 < DOWNLOAD_RETRIES:
                time.sleep(RETRY_DELAY)
            else:
                logging.error("Final download attempt failed.")
                sys.exit(1)

    # --- Verification Phase ---
    logging.info("Verifying file integrity...")
    zip_hash = calculate_sha256(local_zip_filepath)
    logging.info(f"SHA256: {zip_hash}")

    # --- Extraction Phase ---
    # (Same as your original logic)
    temp_extract_dir = data_dir.with_name(data_dir.name + "_tmp")
    try:
        logging.info(f"Extracting data...")
        if temp_extract_dir.exists():
            shutil.rmtree(temp_extract_dir)
        temp_extract_dir.mkdir()

        with zipfile.ZipFile(local_zip_filepath, 'r') as zip_ref:
            file_list = [m.filename for m in zip_ref.infolist() if not m.is_dir()]
            common_path = os.path.commonpath(file_list)

            for member in zip_ref.infolist():
                rel_path = os.path.relpath(member.filename, common_path)
                target_path = temp_extract_dir / rel_path
                if member.is_dir():
                    target_path.mkdir(parents=True, exist_ok=True)
                else:
                    target_path.parent.mkdir(parents=True, exist_ok=True)
                    with zip_ref.open(member) as source, open(target_path, "wb") as target:
                        shutil.copyfileobj(source, target)
        
        if data_dir.exists():
            shutil.rmtree(data_dir)
        shutil.move(str(temp_extract_dir), str(data_dir))
        
        if not args.no_cleanup:
            os.remove(local_zip_filepath)
            logging.info("Cleaned up ZIP file.")

    except Exception as e:
        logging.error(f"Extraction failed: {e}")
        sys.exit(1)

if __name__ == '__main__':
    setup_logging()
    main()