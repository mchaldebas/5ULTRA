import gdown
import zipfile
import argparse
import logging
import hashlib
import shutil
import os
import sys
from pathlib import Path

# --- Configuration ---
GOOGLE_DRIVE_FILE_ID = '1zVYEZt72qgCgye8_Az6zfAmUk1-cGp9U'

def setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        stream=sys.stdout,
    )

def calculate_sha256(filepath: Path) -> str:
    sha256_hash = hashlib.sha256()
    with open(filepath, "rb") as f:
        for byte_block in iter(lambda: f.read(4096), b""):
            sha256_hash.update(byte_block)
    return sha256_hash.hexdigest()

def main():
    parser = argparse.ArgumentParser(description="Download and set up 5ULTRA data.")
    default_data_dir = Path.home() / ".5ULTRA" / "data"
    parser.add_argument('--data-dir', type=Path, default=default_data_dir)
    parser.add_argument('--force', action='store_true')
    parser.add_argument('--no-cleanup', action='store_true')
    args = parser.parse_args()

    data_dir: Path = args.data_dir
    local_zip_filepath = data_dir / 'data.zip'

    # Check if data already exists
    if data_dir.exists() and any(data_dir.iterdir()) and not args.force:
        logging.info(f"Data directory '{data_dir}' already exists. Use --force to overwrite.")
        sys.exit(0)

    data_dir.mkdir(parents=True, exist_ok=True)

    # --- Download Phase ---
    logging.info(f"Downloading data from Google Drive...")
    
    try:
        # gdown handles the confirmation token and progress bar automatically
        gdown.download(
            id=GOOGLE_DRIVE_FILE_ID, 
            output=str(local_zip_filepath), 
            quiet=False, 
            fuzzy=True
        )
    except Exception as e:
        logging.error(f"Download failed: {e}")
        sys.exit(1)

    # --- Verification & Extraction ---
    if not local_zip_filepath.exists() or local_zip_filepath.stat().st_size < 1000000:
        logging.error("Downloaded file is too small. Google Drive may have blocked the request.")
        sys.exit(1)

    logging.info("Extracting data...")
    temp_extract_dir = data_dir.with_name(data_dir.name + "_tmp")
    
    try:
        if temp_extract_dir.exists(): shutil.rmtree(temp_extract_dir)
        temp_extract_dir.mkdir()

        with zipfile.ZipFile(local_zip_filepath, 'r') as zip_ref:
            file_list = [m.filename for m in zip_ref.infolist() if not m.is_dir()]
            common_path = os.path.commonpath(file_list) if file_list else ""

            for member in zip_ref.infolist():
                rel_path = os.path.relpath(member.filename, common_path)
                target_path = temp_extract_dir / rel_path
                if member.is_dir():
                    target_path.mkdir(parents=True, exist_ok=True)
                else:
                    target_path.parent.mkdir(parents=True, exist_ok=True)
                    with zip_ref.open(member) as source, open(target_path, "wb") as target:
                        shutil.copyfileobj(source, target)
        
        if not args.no_cleanup and local_zip_filepath.exists():
            os.remove(local_zip_filepath)

        # 2. Now replace the old directory with the new one
        if data_dir.exists(): 
            shutil.rmtree(data_dir)
        shutil.move(str(temp_extract_dir), str(data_dir))
        
        logging.info("Setup complete.")
            
    except Exception as e:
        logging.error(f"Extraction failed: {e}")
        sys.exit(1)

if __name__ == '__main__':
    setup_logging()
    main()