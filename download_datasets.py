import os
import zipfile
import subprocess

def download_and_extract_zenodo(record_id, zip_filename, extract_folder):
    try:
        print(f"\tLoading from Zenodo repository n°{record_id}...")
        subprocess.run(["zenodo_get", "-g", zip_filename, str(record_id)], check=True)
        print("\tDownload completed.")
        if os.path.exists(zip_filename):
            os.makedirs(extract_folder, exist_ok=True)
            print(f"\tExtraction in progress in the folder : '{extract_folder}'...")
            with zipfile.ZipFile(zip_filename, 'r') as zip_ref:
                zip_ref.extractall(extract_folder)
            print("\tExtraction completed successfully !")
        else:
            print(f"[Error] The file '{zip_filename}' was not found after downloading.")
    except subprocess.CalledProcessError:
        print("[Error] zenodo-get command failed. Please ensure you have the 'zenodo_get' tool installed and properly configured.")
    except Exception as e:
        print(f"[Error] An unexpected error occurred: {e}")

if __name__ == "__main__":
    RECORD_ID = "18965466" 
    ZIP_FILE_NAME = "datasets.zip"
    EXTRACTION_DIR = "datasets"

    download_and_extract_zenodo(RECORD_ID, ZIP_FILE_NAME, EXTRACTION_DIR)