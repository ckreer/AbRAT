#!/usr/bin/env python3

import os
import argparse
import requests
import subprocess
import datetime

# URLs to AURR OGRDB:
URL_IGH_VDJ = "https://ogrdb.airr-community.org/api/germline/set/Human/IGH_VDJ/published/ungapped_ex"
URL_IGKappa_VJ = "https://ogrdb.airr-community.org/api/germline/set/Human/IGKappa_VJ/published/ungapped_ex"
URL_IGLambda_VJ = "https://ogrdb.airr-community.org/api/germline/set/Human/IGLambda_VJ/published/ungapped_ex"

# Lokale Datei-Namen, unter denen wir speichern
FILE_IGH_VDJ = "IGH_VDJ.fasta"
FILE_IGKappa_VJ = "IGKappa_VJ.fasta"
FILE_IGLambda_VJ = "IGLambda_VJ.fasta"

# Output-Files nach Segmenttyp
FILE_HUMAN_V = "human_V.fasta"
FILE_HUMAN_D = "human_D.fasta"
FILE_HUMAN_J = "human_J.fasta"


def download_file(url, local_filename):
    """
    Downloads a file from a URL and saves it locally.

    This function uses the Requests library to download a file in streaming mode from the specified URL
    and writes it to the given local filename.

    Parameters:
        url (str): The URL to download the file from.
        local_filename (str): The local path where the file will be saved.

    Returns:
        None
    """
    print(f"Downloading {url} -> {local_filename}")
    r = requests.get(url, stream=True)
    r.raise_for_status()
    with open(local_filename, 'wb') as f:
        for chunk in r.iter_content(chunk_size=8192):
            f.write(chunk)


def split_fasta_by_segment(input_files, file_v, file_d, file_j):
    """
    Splits FASTA entries from input files into three output files based on segment type.

    This function reads FASTA files provided in 'input_files' and distributes the entries into three separate
    output files (V, D, and J) based on the header content:
      - V: Headers containing "IGHV", "IGKV", or "IGLV".
      - D: Headers containing "IGHD".
      - J: Headers containing "IGHJ", "IGKJ", or "IGLJ".

    Parameters:
        input_files (list): List of file paths to FASTA files.
        file_v (str): Output file path for V segments.
        file_d (str): Output file path for D segments.
        file_j (str): Output file path for J segments.

    Returns:
        None
    """
    # Open V, D, and J files for writing data
    with open(file_v, 'w') as fv, open(file_d, 'w') as fd, open(file_j, 'w') as fj:
        for fasta_file in input_files:
            # Iterate through FASTA files to read data
            with open(fasta_file, 'r') as fin:
                current_entry = []
                current_header = None

                for line in fin:
                    if line.startswith('>'):
                        # If an entry exists, write it first to the correct output file
                        if current_entry and current_header:
                            if is_v_segment(current_header):
                                fv.writelines(current_entry)
                            elif is_d_segment(current_header):
                                fd.writelines(current_entry)
                            elif is_j_segment(current_header):
                                fj.writelines(current_entry)

                        # New sequence: store header and initialize current entry
                        current_header = line.strip()
                        current_entry = [line]  # Header is the first line
                    else:
                        # Append sequence line
                        current_entry.append(line)

                # Write the last entry in the file, if any
                if current_entry and current_header:
                    if is_v_segment(current_header):
                        fv.writelines(current_entry)
                    elif is_d_segment(current_header):
                        fd.writelines(current_entry)
                    elif is_j_segment(current_header):
                        fj.writelines(current_entry)


def is_v_segment(header):
    """
    Determines if a FASTA header corresponds to a V segment.

    This function checks if the header contains any of the substrings: "IGHV", "IGKV", or "IGLV".

    Parameters:
        header (str): The FASTA header string.

    Returns:
        bool: True if the header indicates a V segment, otherwise False.
    """
    return ("IGHV" in header) or ("IGKV" in header) or ("IGLV" in header)


def is_d_segment(header):
    """
    Determines if a FASTA header corresponds to a D segment.

    This function checks if the header contains the substring "IGHD".

    Parameters:
        header (str): The FASTA header string.

    Returns:
        bool: True if the header indicates a D segment, otherwise False.
    """
    return "IGHD" in header


def is_j_segment(header):
    """
    Determines if a FASTA header corresponds to a J segment.

    This function checks if the header contains any of the substrings: "IGHJ", "IGKJ", or "IGLJ".

    Parameters:
        header (str): The FASTA header string.

    Returns:
        bool: True if the header indicates a J segment, otherwise False.
    """
    return ("IGHJ" in header) or ("IGKJ" in header) or ("IGLJ" in header)


def make_blast_db(input_fasta, output_name):
    """
    Creates a nucleotide BLAST database from an input FASTA file using makeblastdb.

    This function runs the makeblastdb command with the appropriate parameters to generate a nucleotide
    database. The sequence identifiers are parsed from the FASTA file.

    Parameters:
        input_fasta (str or Path): Path to the input FASTA file.
        output_name (str): Name of the output BLAST database.

    Returns:
        None
    """
    print(f"Creating BLAST DB: {output_name}")
    subprocess.run([
        "makeblastdb",
        "-parse_seqids",
        "-dbtype", "nucl",
        "-in", input_fasta,
        "-out", output_name
    ], check=True)


def main():
    parser = argparse.ArgumentParser(description="Sequence annotation and quality check")
    parser.add_argument('--output_folder', type=str, required=True, help='Output Folder')

    args = parser.parse_args()

    # 1) Check output folder
    if not os.path.exists(args.output_folder):
        os.makedirs(args.output_folder)

    # 2) Download AIRR data
    download_file(URL_IGH_VDJ, os.path.join(args.output_folder, FILE_IGH_VDJ))
    download_file(URL_IGKappa_VJ, os.path.join(args.output_folder, FILE_IGKappa_VJ))
    download_file(URL_IGLambda_VJ, os.path.join(args.output_folder, FILE_IGLambda_VJ))

    # 3) Split by V, D, J
    split_fasta_by_segment(
        [os.path.join(args.output_folder, FILE_IGH_VDJ),
         os.path.join(args.output_folder, FILE_IGKappa_VJ),
         os.path.join(args.output_folder, FILE_IGLambda_VJ)],
        os.path.join(args.output_folder, FILE_HUMAN_V),
        os.path.join(args.output_folder, FILE_HUMAN_D),
        os.path.join(args.output_folder, FILE_HUMAN_J)
    )

    # 4) makeblastdb
    make_blast_db(os.path.join(args.output_folder, FILE_HUMAN_V), os.path.join(args.output_folder, "human_V"))
    make_blast_db(os.path.join(args.output_folder, FILE_HUMAN_D), os.path.join(args.output_folder, "human_D"))
    make_blast_db(os.path.join(args.output_folder, FILE_HUMAN_J), os.path.join(args.output_folder, "human_J"))

    # 5) Save a time stamp
    with open(os.path.join(args.output_folder, "TIME_STAMP.txt"), 'w') as ts_file:
        ts_file.write(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        ts_file.write("\nDatabase created with makedb and AIRR data at:\n")
        ts_file.write("\n".join([URL_IGH_VDJ, URL_IGKappa_VJ, URL_IGLambda_VJ]))

    print("All done. Created human_V, human_D, human_J and corresponding BLAST DBs.")

if __name__ == "__main__":
    main()
