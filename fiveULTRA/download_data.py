import requests
import os
import zipfile
import argparse

def main():
    parser = argparse.ArgumentParser(
        description="Download necessary data for 5ULTRA"
    )
    parser.add_argument('--data-dir', type=str, default=os.path.join(os.path.expanduser("~"), ".5ULTRA", "data"), help='Path to the data directory')
    args = parser.parse_args()
    data_dir = args.data_dir
    url = 'https://www.dropbox.com/scl/fo/pj5ls1m8f9gvc6b13tqle/AHHge_c_7gZHtn15ogomSOQ?rlkey=n233zz0nbmz1z9qhyihspkcmi&st=j7xz1038&dl=1'
    local_filename = 'data.zip'

    # Ensure the target directory exists
    os.makedirs(data_dir, exist_ok=True)
    local_filepath = os.path.join(data_dir, local_filename)

    # Download the file
    print(f"Downloading data to {local_filepath}...")
    with requests.get(url, stream=True) as r:
        r.raise_for_status()
        with open(local_filepath, 'wb') as f:
            for chunk in r.iter_content(chunk_size=8192):
                f.write(chunk)
    print("Download complete.")

    # Extract the ZIP file
    print("Extracting data...")
    try:
        with zipfile.ZipFile(local_filepath, 'r') as zip_ref:
            zip_ref.extractall(data_dir)
        print(f"Data extracted to {data_dir}")
    except zipfile.BadZipFile:
        print("Error: The downloaded file is not a valid ZIP archive.")
    
    # Clean up the ZIP file after extraction
    os.remove(local_filepath)
    print("Cleanup complete. Data is ready.")

if __name__ == '__main__':
    main()
