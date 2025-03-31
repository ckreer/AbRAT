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
    """Einfaches Download einer Datei via Requests."""
    print(f"Downloading {url} -> {local_filename}")
    r = requests.get(url, stream=True)
    r.raise_for_status()
    with open(local_filename, 'wb') as f:
        for chunk in r.iter_content(chunk_size=8192):
            f.write(chunk)


def split_fasta_by_segment(input_files, file_v, file_d, file_j):
    """
    Liest die FASTA-Dateien in 'input_files' ein und verteilt die Einträge
    in drei Ausgabedateien (V, D, J).
    Wir erkennen den Segmenttyp an den Headern:
    - V: >IGHV, >IGKV, >IGLV
    - D: >IGHD
    - J: >IGHJ, >IGKJ, >IGLJ
    """
    # open v, d and j file to write data
    with open(file_v, 'w') as fv, open(file_d, 'w') as fd, open(file_j, 'w') as fj:
        for fasta_file in input_files:
            # go through fasta files to read data
            with open(fasta_file, 'r') as fin:
                current_entry = []
                current_header = None

                for line in fin:
                    if line.startswith('>'):
                        # if any entry, write it first
                        if current_entry and current_header:
                            if is_v_segment(current_header):
                                fv.writelines(current_entry)
                            elif is_d_segment(current_header):
                                fd.writelines(current_entry)
                            elif is_j_segment(current_header):
                                fj.writelines(current_entry)

                        # New seqeunce, remember header
                        current_header = line.strip()
                        current_entry = [line]  # header aas first line
                    else:
                        # sequence line
                        current_entry.append(line)

                # last entry in file?
                if current_entry and current_header:
                    if is_v_segment(current_header):
                        fv.writelines(current_entry)
                    elif is_d_segment(current_header):
                        fd.writelines(current_entry)
                    elif is_j_segment(current_header):
                        fj.writelines(current_entry)


def is_v_segment(header):
    # V-Segmente erkennen an "IGHV", "IGKV", "IGLV"
    # (header könnte z.B. so aussehen: ">IGHV1-2*01 Some description")
    return ("IGHV" in header) or ("IGKV" in header) or ("IGLV" in header)


def is_d_segment(header):
    # D-Segmente nur "IGHD"
    return ("IGHD" in header)


def is_j_segment(header):
    # J-Segmente: "IGHJ", "IGKJ", "IGLJ"
    return ("IGHJ" in header) or ("IGKJ" in header) or ("IGLJ" in header)


def make_blast_db(input_fasta, output_name):
    """ Runs makeblastdb to generate a nukleotide database. """
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
