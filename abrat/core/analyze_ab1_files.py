#!/usr/bin/env python3

import argparse
import subprocess
import os
import logging
import pandas as pd
import numpy as np

from pathlib import Path
from Bio import SeqIO

from abrat.core.utils import (decompress,
                              get_indels,
                              get_sequence,
                              translate_nt_to_aa,
                              generate_q_check_df,
                              get_isotype_with_blast,
                              make_quality_check_excel,
                              log_message, date_stamp,
                              str2bool)

# Settings
path_to_human_V = "/app/data/database/igblastdb/human_V"
path_to_human_D = "/app/data/database/igblastdb/human_D"
path_to_human_J = "/app/data/database/igblastdb/human_J"

airr_fmt = "19"

def combine_ab1_files(input_path, output_path, cf_name):
    """
    Combines .ab1 files into a single FASTA file.

    This function reads all *.ab1 files from the specified input directory,
    extracts the record names and sequences using SeqIO, concatenates them into
    a single FASTA-formatted string, and writes the result to an output file.

    Parameters:
        input_path (str or Path): Path to the directory containing .ab1 files.
        output_path (str or Path): Path to the directory where the output FASTA file will be saved.
        cf_name (str): Name of the output FASTA file.

    Returns:
        str: Full path to the combined FASTA file if the operation is successful.
        tuple: (False, error_message) if an error occurs.
    """
    try:
        # Create output directory if it does not exist
        if not os.path.exists(output_path):
            os.makedirs(output_path)

        combined = ""

        ab1_files = list(Path(input_path).rglob("*.ab1"))

        if not ab1_files:
            raise FileNotFoundError("No .ab1-files found in " + str(input_path))

        for ab1_file in ab1_files:
            with ab1_file.open("rb") as input_file:
                # Append record name and sequence to the combined string in FASTA format
                for record in SeqIO.parse(input_file, 'abi'):
                    combined += ">" + record.name + "\n" + str(record.seq) + "\n"

        with open(os.path.join(output_path, cf_name), 'w') as c_fasta:
            c_fasta.write(combined)

        return os.path.join(output_path, cf_name)

    except FileNotFoundError as e:
        return False, f"Error: {e}"
    except Exception as e:
        return False, f"General error: {e}"


def run_igblast(path_to_input_file, output_folder, outformat):
    """
    Runs IgBLAST as a subprocess using the specified input file and captures its stdout and stderr.

    This function verifies the existence of the input file and output folder, constructs the IgBLAST command,
    and executes it using the subprocess module. If IgBLAST runs successfully, the function returns True
    along with the stdout output. In case of an error, it returns False and the corresponding error message.

    Parameters:
        path_to_input_file (str or Path): Path to the input file.
        output_folder (str or Path): Path to the output folder.
        outformat (str): Output format option for IgBLAST.

    Returns:
        tuple: A tuple (bool, str) where the first element is True if IgBLAST ran successfully, otherwise False.
               The second element is the stdout output on success, or an error message on failure.
    """
    logging.basicConfig(level=logging.INFO,
                        format='%(asctime)s - %(levelname)s - %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')
    logger = logging.getLogger()

    if not os.path.exists(path_to_input_file):
        logger.error(f"Input folder does not exist: {path_to_input_file}")
        return False, f"Input folder does not exist: {path_to_input_file}"

    if not os.path.exists(output_folder):
        os.makedirs(output_folder)

    # Construct the IgBLAST command
    cmd = [
        "igblastn",
        "-germline_db_V", path_to_human_V,
        "-germline_db_D", path_to_human_D,
        "-germline_db_J", path_to_human_J,
        "-auxiliary_data", "optional_file/human_gl.aux",
        "-query", path_to_input_file,
        "-outfmt", outformat,
        "-extend_align5end"
    ]
    try:
        print(log_message("Running IgBLAST"), flush=True)
        result = subprocess.run(cmd, check=True, capture_output=True, text=True)
        return True, result.stdout

    except subprocess.CalledProcessError as e:
        msg = f"Error from IgBLAST: {e}\nStderr:\n{e.stderr}"
        logger.error(msg)
        return False, msg

    except Exception as e:
        msg = f"General error: {e}"
        logger.error(msg)
        return False, msg

def get_infos_from_igblast(blast_result):
    """
    Extracts information from an IgBLAST result string.

    This function parses the IgBLAST result string to extract the sample name and various features
    such as gene segments, sequence regions, and alignment details. Note that FWR4_START is defined
    relative to the J gene while FWR4_END (like all other positions) is relative to the whole sequence;
    however, FWR4_START is exported with reference to the original sequence.

    Parameters:
        blast_result (str): The IgBLAST result as a string.

    Returns:
        tuple: A tuple (name, features) where 'name' is the sample name (str) and 'features' is a dictionary
               containing extracted features.
    """
    # Workaround for getting the correct C-terminus
    c_term_offset = 0
    fwr4_correction = 0

    parts = blast_result.split('# ')
    # Initialize values for quality check
    name, orientation, fwr1_start, fwr1_found, j_start, j_end, j_found, v_length, q_overhang = \
        "N/A", "N/A", -1, False, -1, -1, False, 0, 0
    # q_overhang: if the query does not start at 1, add the unknown bases before
    # Initialize values for other important BLAST results
    v_gene, j_gene, d_gene, v_btop, j_btop = "N/A", "N/A", "N/A", "N/A", "N/A"
    chain_type, stop_codon, frame, productive = "N/A", "N/A", "N/A", "N/A"
    cdr3_nt, cdr3_aa, cdr3_start, cdr3_end = "N/A", "N/A", "N/A", "N/A"
    j_nt, fwr4_start, fwr4_end, fwr4_nt, fwr4_aa, fwr4_frame, fwr4_stop_codon = \
        "N/A", "N/A", 0, "N/A", "N/A", "N/A", "N/A"
    fwr1_end, cdr1_start, cdr1_end, fwr2_start, fwr2_end = "N/A", "N/A", "N/A", "N/A", "N/A"
    cdr2_start, cdr2_end, fwr3_start, fwr3_end, v_identity = "N/A", "N/A", "N/A", "N/A", "N/A"
    cdr1_found, fwr2_found, cdr2_found, fwr3_found, cdr3_found, fwr4_found = False, False, False, False, False, False
    truncated_kappa = False

    for part in parts:
        if part.startswith('Query: '):
            name = part.split(': ')[1].rstrip('\n')
        elif part.startswith('V-(D)-J rearrangement'):
            # Light chain specific:
            if ("VK\t" in part.split('\n')[1]) or ("VL\t" in part.split('\n')[1]):
                v_gene = part.split('\n')[1].split('\t')[0]
                j_gene = part.split('\n')[1].split('\t')[1]
                chain_type = part.split('\n')[1].split('\t')[2]
                stop_codon = part.split('\n')[1].split('\t')[3]
                frame = part.split('\n')[1].split('\t')[4]
                productive = part.split('\n')[1].split('\t')[5]
            # Heavy chain specific:
            if "VH\t" in part.split('\n')[1]:
                v_gene = part.split('\n')[1].split('\t')[0]
                d_gene = part.split('\n')[1].split('\t')[1]
                j_gene = part.split('\n')[1].split('\t')[2]
                chain_type = part.split('\n')[1].split('\t')[3]
                stop_codon = part.split('\n')[1].split('\t')[4]
                frame = part.split('\n')[1].split('\t')[5]
                productive = part.split('\n')[1].split('\t')[6]
            # For both chain types, the second last part is orientation and the last is frame shift
            orientation = part.split('\n')[1].split('\t')[-2]
        elif part.startswith('Sub-region sequence'):
            cdr3_nt = part.split('\n')[1].split('\t')[1]
            cdr3_aa = part.split('\n')[1].split('\t')[2]
            cdr3_start = int(part.split('\n')[1].split('\t')[3])
            cdr3_end = int(part.split('\n')[1].split('\t')[4])
            cdr3_found = True
        elif part.startswith('Alignment summary'):
            subparts = part.split('\n')
            for subpart in subparts:
                if subpart.startswith('FR1'):
                    fwr1_start = int(subpart.split('\t')[1])
                    fwr1_end = int(subpart.split('\t')[2])
                    fwr1_found = True
                if subpart.startswith('CDR1'):
                    cdr1_start = int(subpart.split('\t')[1])
                    cdr1_end = int(subpart.split('\t')[2])
                    cdr1_found = True
                if subpart.startswith('FR2'):
                    fwr2_start = int(subpart.split('\t')[1])
                    fwr2_end = int(subpart.split('\t')[2])
                    fwr2_found = True
                if subpart.startswith('CDR2'):
                    cdr2_start = int(subpart.split('\t')[1])
                    cdr2_end = int(subpart.split('\t')[2])
                    cdr2_found = True
                if subpart.startswith('FR3'):
                    fwr3_start = int(subpart.split('\t')[1])
                    fwr3_end = int(subpart.split('\t')[2])
                    fwr3_found = True
                elif subpart.startswith('Total'):
                    v_length = int(subpart.split('\t')[3])
                    v_identity = float(subpart.split('\t')[7])
        elif 'hits found' in part and not part.startswith("0"):
            v_start = int(part.split('\n')[1].split('\t')[10])  # start alignment of V gene (should be 1 if complete match)
            s_start = int(part.split('\n')[1].split('\t')[8])   # first nucleotide of sequence matching to V gene
            q_overhang = v_start - 1 if s_start > v_start else 0  # add overhang if original sequence is longer than V start
            v_btop = part.split('\n')[1].split('\t')[16]
            # Check if FR1 is complete
            if fwr1_found and (v_start > 1) and (s_start < v_start):
                if (v_start - s_start < 20) and (chain_type == "VK"):
                    fwr1_found = True
                    truncated_kappa = True
                else:
                    fwr1_found = False
            for subpart in part.split('\n'):
                if subpart.startswith('J') and not j_found:
                    j_start = int(subpart.split('\t')[8])
                    j_end = int(subpart.split('\t')[9])
                    j_nt = subpart.split('\t')[14]
                    j_btop = subpart.split('\t')[16]
                    if cdr3_found and (cdr3_end != "N/A"):
                        if j_start > cdr3_end:  # modified to account for insertions at the end of CDR3
                            fwr4_start = 0
                        else:
                            fwr4_start = cdr3_end - j_start + 1
                        fwr4_end = j_end
                        if (fwr4_start >= 0) and (fwr4_end > 0) and (j_btop != "N/A"):
                            fwr4_nt = j_nt[fwr4_start:]
                            if len(fwr4_nt) > 0:
                                # Workaround to solve IgBLAST problem where C/G is sometimes missed or not at the 5' end
                                if len([fwr4_nt[i:i+3] for i in range(0, len(fwr4_nt), 3)][-1]) == 1:
                                    fwr4_correction = -1
                                    c_term_offset = -1
                                # If the last triplet is missing one nucleotide, add 1 nucleotide from the original sequence
                                elif len([fwr4_nt[i:i+3] for i in range(0, len(fwr4_nt), 3)][-1]) == 2:
                                    fwr4_correction = +1
                                    c_term_offset = +1

                                # Perform frame analysis:
                                ins, dels = get_indels(decompress(j_btop)[fwr4_start:], decompression=False)
                                if ((ins - dels) % 3) == 0:
                                    fwr4_frame = "In-frame"
                                else:
                                    fwr4_frame = "Out-of-frame"

                                fwr4_aa = translate_nt_to_aa(fwr4_nt)
                                fwr4_stop_codon = "Yes" if "*" in fwr4_aa else "No"
                                fwr4_found = True

                    j_found = True
        # Elongate 5' end if BLAST did not start alignment at 1 but more sequence is available
        if fwr1_start > q_overhang:
            v_length += q_overhang
            fwr1_start -= q_overhang

    if fwr1_found and j_found:
        total_length = j_end - fwr1_start if j_end > fwr1_start else 0

    orient_found = False if orientation == "N/A" else True
    v_length_found = False if v_length < 1 else True

    # Check if the V alignment is complete at the 5' end
    if (q_overhang > 0) or truncated_kappa or (fwr1_start == -1):
        full_v_alignment = False
    else:
        full_v_alignment = True

    # Add constant region here
    # TODO: IgBLAST sometimes adds the G to the J gene, sometimes to the C region.
    # Determine when this occurs to obtain the correct start (this is a workaround).
    if fwr4_found:
        constant_region_start = fwr4_end + 1 + c_term_offset
        fwr4_start = fwr4_start + j_start
    else:
        constant_region_start = False
    igblast_pass = fwr1_found and orient_found and v_length_found and j_found and fwr3_found

    features = {
        'V_GENE': v_gene,
        'TOP_V': v_gene.split('*')[0],
        'D_GENE': d_gene,
        'TOP_D': d_gene.split('*')[0],
        'J_GENE': j_gene,
        'TOP_J': j_gene.split('*')[0],
        'V_IDENTITY': v_identity,
        'CHAIN_TYPE': chain_type,
        'STOP_CODON': stop_codon,
        'FRAME': frame,
        'PRODUCTIVE': productive,
        'ORIENTATION': orientation,
        'FWR1_START': fwr1_start,
        'FWR1_END': fwr1_end,
        'CDR1_START': cdr1_start,
        'CDR1_END': cdr1_end,
        'FWR2_START': fwr2_start,
        'FWR2_END': fwr2_end,
        'CDR2_START': cdr2_start,
        'CDR2_END': cdr2_end,
        'FWR3_START': fwr3_start,
        'FWR3_END': fwr3_end,
        'CDR3_START': cdr3_start,
        'CDR3_END': cdr3_end,
        'CDR3_NT': cdr3_nt,
        'CDR3_AA': cdr3_aa,
        'V_BTOP': v_btop,
        'J_START': j_start,
        'J_END': j_end,
        'J_BTOP': j_btop,
        'FWR4_START': fwr4_start,
        'FWR4_END': fwr4_end,
        'FWR4_CORRECTION': fwr4_correction,
        'FWR4_NT': fwr4_nt,
        'FWR4_AA': fwr4_aa,
        'FWR4_FRAME': fwr4_frame,
        'FWR4_STOP_CODON': fwr4_stop_codon,
        'C_START': constant_region_start,
        'V_LENGTH': v_length,
        'FWR1_FOUND': fwr1_found,
        'J_FOUND': j_found,
        'ORIENT_FOUND': orient_found,
        'LENGTH_FOUND': v_length_found,
        'FULL_V_ALIGNMENT': full_v_alignment,
        'IGBLAST_PASSED': igblast_pass
    }
    return name, features


def combine_ig_blast_and_q_check_data(q_check_df, igblast_dict, combi_ex_name, out_dir, isotype_determination):
    """
    Combines IgBLAST data with quality check data into a single DataFrame and exports the result as an Excel file.

    This function merges a quality check DataFrame with IgBLAST results based on sample names,
    computes sequence segments from the original sequence, optionally determines isotypes,
    and generates an output Excel file containing the combined data.

    Parameters:
        q_check_df (pd.DataFrame): DataFrame containing quality check data.
        igblast_dict (dict): Dictionary containing IgBLAST results keyed by sample name.
        combi_ex_name (str): Name of the output Excel file.
        out_dir (str or Path): Directory where the output Excel file will be saved.
        isotype_determination (bool): Flag indicating whether isotype determination should be performed.

    Returns:
        pd.DataFrame: Combined DataFrame with merged and computed information.
    """
    # Create a DataFrame from igblast_dict, using the dictionary keys as sample names
    igblast_df = pd.DataFrame.from_dict(igblast_dict, orient='index').reset_index()
    igblast_df.rename(columns={'index': 'SAMPLE_NAME'}, inplace=True)

    # Add a SAMPLE_NAME column to q_check_df by removing the '.ab1' suffix from FILE_NAME
    q_check_df['SAMPLE_NAME'] = q_check_df['FILE_NAME'].str.replace('.ab1', '')

    # Select a subset of columns from q_check_df
    subset_q_check = q_check_df[['SAMPLE_NAME', 'SUBSET', 'TISSUE', 'CHAIN_PCR',
                                 'MEAN_PHRED_PASSED', 'INNER_N_PASSED',
                                 'LENGTH_PASSED', 'INNER_N', 'QCHECK_PASSED', 'ORIGINAL_SEQ',
                                 'TRIMMED_SEQ', 'MASKED_SEQ', 'COHORT', 'SUBJECT',
                                 'TIME_POINT', 'PLATE', 'PRIMER_SET', 'WELL',
                                 'SOURCE', 'SUBSOURCE',
                                 # Add any additional columns as needed
                                 ]]

    # Merge the quality check data with the IgBLAST DataFrame on SAMPLE_NAME
    combined_df = pd.merge(subset_q_check, igblast_df, on='SAMPLE_NAME', how='left')

    # Calculate sequence segments based on the original sequence.
    # Assumption: get_sequence and translate_nt_to_aa are imported and functioning.
    combined_df['FWR1_NT'] = combined_df.apply(
        lambda row: get_sequence(row['FWR1_START'], row['FWR1_END'], row['ORIGINAL_SEQ'], -1), axis=1)
    combined_df['CDR1_NT'] = combined_df.apply(
        lambda row: get_sequence(row['CDR1_START'], row['CDR1_END'], row['ORIGINAL_SEQ'], -1), axis=1)
    combined_df['CDR1_AA'] = combined_df['CDR1_NT'].apply(translate_nt_to_aa)
    combined_df['FWR2_NT'] = combined_df.apply(
        lambda row: get_sequence(row['FWR2_START'], row['FWR2_END'], row['ORIGINAL_SEQ'], -1), axis=1)
    combined_df['CDR2_NT'] = combined_df.apply(
        lambda row: get_sequence(row['CDR2_START'], row['CDR2_END'], row['ORIGINAL_SEQ'], -1), axis=1)
    combined_df['CDR2_AA'] = combined_df['CDR2_NT'].apply(translate_nt_to_aa)
    combined_df['FWR3_NT'] = combined_df.apply(
        lambda row: get_sequence(row['FWR3_START'], row['FWR3_END'], row['ORIGINAL_SEQ'], -1), axis=1)

    # FWR4: If a correction is required, calculate FWR4_NT and FWR4_AA
    if combined_df['FWR4_CORRECTION'].iloc[0]:
        combined_df['FWR4_NT'] = combined_df.apply(
            lambda row: get_sequence(row['FWR4_START'], int(row['FWR4_END']) + int(row['FWR4_CORRECTION']),
                                     row['ORIGINAL_SEQ'], -1),
            axis=1
        )
        combined_df['FWR4_AA'] = combined_df['FWR4_NT'].apply(translate_nt_to_aa)
    else:
        combined_df['FWR4_NT'] = None
        combined_df['FWR4_AA'] = None

    # Compute length information for CDR3
    combined_df['CDR3_NT_LENGTH'] = combined_df['CDR3_NT'].apply(lambda x: len(x) if x != "N/A" else "N/A")
    combined_df['CDR3_AA_LENGTH'] = combined_df['CDR3_AA'].apply(lambda x: len(x) if x != "N/A" else "N/A")

    # Calculate C_NT and C_AA (assumes that C_START is available in the IgBLAST result)
    combined_df['C_NT'] = combined_df.apply(
        lambda row: get_sequence(row['C_START'], len(row['ORIGINAL_SEQ']) - 20, row['ORIGINAL_SEQ'], -1), axis=1)
    combined_df['C_AA'] = combined_df['C_NT'].apply(translate_nt_to_aa)

    # Isotype determination
    if isotype_determination:
        print(log_message("Determining Isotypes"), flush=True)
        combined_df['ISOTYPE'] = get_isotype_with_blast(combined_df['C_NT'], out_dir)
        combined_df['TOP_ISOTYPE'] = combined_df['ISOTYPE'].apply(lambda x: x.split("*")[0])
    else:
        combined_df['ISOTYPE'] = "N.D."
        combined_df['TOP_ISOTYPE'] = "N.D."

    # If desired: Correct the isotype column for light chains
    change_index = combined_df[combined_df['CHAIN_PCR'] != "HC"].index
    combined_df.loc[change_index, "ISOTYPE"] = np.nan

    # Generate a B_CELL_ID using several columns
    combined_df['B_CELL_ID'] = combined_df[
        ['COHORT', 'SUBJECT', 'TIME_POINT', 'TISSUE', 'SUBSET', 'PLATE', 'WELL']].apply(
        lambda x: "_".join(x.astype(str)), axis=1
    )

    # Reorder columns (e.g., 'SELECT_HC', 'SELECT_KC', 'SELECT_LC', etc.)
    features = ['B_CELL_ID', 'SAMPLE_NAME', 'COHORT', 'SUBJECT', 'TIME_POINT',
                'TISSUE', 'SUBSET', 'PLATE', 'WELL', 'PRIMER_SET', 'SOURCE', 'SUBSOURCE', 'IGBLAST_PASSED', 'V_GENE',
                'TOP_V', 'D_GENE', 'TOP_D', 'J_GENE', 'TOP_J', 'V_IDENTITY', 'CHAIN_PCR', 'STOP_CODON', 'FRAME',
                'PRODUCTIVE', 'ORIENTATION', 'FULL_V_ALIGNMENT', 'FWR1_START', 'FWR1_END', 'FWR1_NT', 'CDR1_START',
                'CDR1_END', 'CDR1_NT', 'FWR2_START', 'FWR2_END', 'FWR2_NT', 'CDR2_START', 'CDR2_END', 'CDR2_NT',
                'FWR3_START', 'FWR3_END', 'FWR3_NT', 'CDR3_NT', 'CDR3_NT_LENGTH', 'CDR3_AA', 'CDR3_AA_LENGTH', 'V_BTOP',
                'J_START', 'J_END', 'J_BTOP', 'FWR4_START', 'FWR4_END', 'FWR4_CORRECTION', 'FWR4_NT', 'FWR4_AA',
                'FWR4_FRAME', 'FWR4_STOP_CODON', 'C_START', 'C_NT', 'C_AA', 'ISOTYPE', 'TOP_ISOTYPE', 'V_LENGTH',
                'LENGTH_PASSED', 'INNER_N', 'INNER_N_PASSED', 'MEAN_PHRED_PASSED', 'QCHECK_PASSED', 'TRIMMED_SEQ',
                'MASKED_SEQ', 'ORIGINAL_SEQ']
    combined_df = combined_df.reindex(columns=features)

    # Export the combined DataFrame as an Excel file
    path_to_combined_excel = Path(out_dir).joinpath(combi_ex_name)
    with pd.ExcelWriter(path_to_combined_excel) as writer:
        combined_df.to_excel(writer, sheet_name="Full_Information", index=False)

    print(log_message(f"{combi_ex_name} exported to {out_dir}"), flush=True)
    return combined_df

def main():
    parser = argparse.ArgumentParser(description="Sequence annotation and quality check")
    parser.add_argument('--input', type=str, required=True, help='Path to input files')
    parser.add_argument('--output', type=str, required=True, help='Path to output')
    parser.add_argument('--project_name', type=str, default="NEW-PROJECT", help='Project name')
    parser.add_argument('--q_cut_off', type=int, default=16, help='Minimum Phred score for base calls')
    parser.add_argument('--mean_q_cut_off', type=int, default=28, help='Minimum mean Phred score per sequence')
    parser.add_argument('--min_length', type=int, default=240, help='Minimum sequence length')
    parser.add_argument('--inner_n', type=int, default=15, help='Maximum number of uncertain bases')
    parser.add_argument('--isotype_determination', type=str2bool, default=False, help='Determine chain isotype')
    parser.add_argument('--airr_export', type=str2bool, default=False, help='Export a copy in airr format')
    args = parser.parse_args()

    print(log_message("Combining all ab1-files from"), args.input, flush=True)
    # combine all ab1 files
    combined_fasta_file = date_stamp(args.project_name+"_combined-sequences.fasta")
    path_to_combined_fasta = combine_ab1_files(args.input, args.output, combined_fasta_file)

    if args.airr_export:
        airr_success, airr_message = run_igblast(path_to_combined_fasta, args.output, airr_fmt)

        if airr_success:
            print(log_message("Generating AIRR output"), flush=True)
            path_to_airr_output_file = os.path.join(args.output,
                                                       date_stamp(args.project_name + "_igblast-annotation-airr.tsv"))
            with open(path_to_airr_output_file, "w") as airr_output_file:
                airr_output_file.write(airr_message)
        else:
            print(log_message("AIRR annotation failed"), flush=True)


    # run igblast for qcheck
    success, message = run_igblast(path_to_combined_fasta, args.output, "7 std qseq sseq btop")

    # save output to file
    path_to_igblast_output_file = os.path.join(args.output, date_stamp(args.project_name+"_igblast-annotation.fmt7"))
    with open(path_to_igblast_output_file, "w") as igblast_output_file:
        igblast_output_file.write(message)

    if success:
        # continue only, if igblast worked
        igblast_results = message.split('# IGBLASTN')[1:]

        igblast_dictionary = {}
        for entry in igblast_results:
            key, value = get_infos_from_igblast(entry)
            igblast_dictionary[key] = value

        # Check ab1 files, igblast infos and save in quality check dataframe
        print(log_message("Performing quality check"), flush=True)
        quality_check_df = generate_q_check_df(args.input,
                                               args.q_cut_off,
                                               args.mean_q_cut_off,
                                               args.min_length,
                                               args.inner_n,
                                               igblast_dictionary)
        # Generate Excel quality check output
        qcheck_excel_name = date_stamp(args.project_name+"_qc-summary.xlsx")
        path_to_qc_excel_file = Path(args.output).joinpath(qcheck_excel_name)

        # Export excel
        qc_excel_success, qc_excel_msg = make_quality_check_excel(quality_check_df, path_to_qc_excel_file)
        if not qc_excel_success:
            print(log_message(qc_excel_msg), flush=True)
            exit(1)
        print(log_message(qc_excel_msg), flush=True)

        # make a composite dataframe from quality_check_df and igblast dictionary
        print(log_message('Preparing summary output with all data'), flush=True)
        combined_excel_name = date_stamp(args.project_name+"_all-sequences.xlsx")
        composite_df = combine_ig_blast_and_q_check_data(quality_check_df,
                                                         igblast_dictionary,
                                                         combined_excel_name,
                                                         args.output,
                                                         args.isotype_determination)

        exit(0)
    else:
        print(log_message("IgBLAST failed."), flush=True)
        exit(1)

if __name__ == "__main__":
    main()
