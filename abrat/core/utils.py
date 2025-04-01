import itertools
import argparse
from subprocess import Popen, PIPE
import os
import sys
import numpy as np
import pandas as pd
import datetime
import matplotlib.pyplot as plt
from pathlib import Path
import math

from Bio import SeqIO

import altair as alt
import plotly.io as pio


def apply_theme_configurations():
    """
    Applies global theme configurations for Matplotlib, Altair, and Plotly.

    This function sets default styling parameters for Matplotlib, defines and enables a custom
    Altair theme, and creates a custom Plotly template that is then set as the default.

    Returns:
        None
    """
    # Matplotlib settings
    plt.rcParams['font.family'] = 'DejaVu Sans'
    plt.rcParams['font.size'] = 12
    plt.rcParams['lines.linewidth'] = 0.5
    plt.rcParams['axes.linewidth'] = 0.5
    plt.rcParams['xtick.major.width'] = 0.5
    plt.rcParams['ytick.major.width'] = 0.5
    plt.rcParams['xtick.minor.width'] = 0.5
    plt.rcParams['ytick.minor.width'] = 0.5

    # Altair: define and enable the custom theme
    def my_altair_theme():
        return {
            "config": {
                "title": {
                    "font": "DejaVu Sans",
                    "fontSize": 12,
                    "fontColor": "black"
                },
                "axis": {
                    "labelFont": "DejaVu Sans",
                    "labelFontSize": 12,
                    "titleFont": "DejaVu Sans",
                    "titleFontSize": 12,
                    "domainColor": "black",
                    "domainWidth": 0.5,
                    "tickColor": "black",
                    "tickSize": 5
                },
                "legend": {
                    "labelFont": "DejaVu Sans",
                    "labelFontSize": 12,
                    "titleFont": "DejaVu Sans",
                    "titleFontSize": 12,
                    "symbolSize": 100
                }
            }
        }

    alt.themes.register("my_altair_theme", my_altair_theme)
    alt.themes.enable("my_altair_theme")

    # Plotly: define and set a custom template as default
    pio.templates["my_plotly_template"] = pio.templates["plotly_white"]
    pio.templates["my_plotly_template"].layout.font.family = "DejaVu Sans"
    pio.templates["my_plotly_template"].layout.font.size = 12
    pio.templates["my_plotly_template"].layout.xaxis.linewidth = 0.5
    pio.templates["my_plotly_template"].layout.yaxis.linewidth = 0.5
    pio.templates.default = "my_plotly_template"

# global settings
# Colors
green = "#caebca"
red = "#FA8072"
orange = "#ffa95b"
light_grey = "#D3D3D3"

# General Functions #
def str2bool(v):
    """
    Converts a string representation of a boolean to a boolean value.

    Parameters:
        v (str or bool): The value to convert.

    Returns:
        bool: True if the value represents a true value; False if it represents a false value.

    Raises:
        argparse.ArgumentTypeError: If the input cannot be interpreted as a boolean.
    """
    if isinstance(v, bool):
        return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')


def timestamp_filename(original_filename, insert_text):
    """
    Generates a new filename by inserting a timestamp between the base name and the extension.
    The timestamp format is _clustered-YYMMDD-HHMMSS.

    Parameters:
        original_filename (str): The original filename (e.g., "myfile.xlsx").
        insert_text (str): Text to insert before the timestamp.

    Returns:
        str: New filename with the timestamp inserted
             (e.g., "myfile_insert_text-230308-153045.xlsx").
    """
    base, ext = os.path.splitext(original_filename)
    timestamp = datetime.datetime.now().strftime("-%y%m%d-%H%M%S")
    return f"{base}_{insert_text}{timestamp}{ext}"


def date_stamp(s):
    """
    Appends the current date to a given string.

    Parameters:
        s (str): The string to which the date will be appended.

    Returns:
        str: A string with the current ISO date followed by an underscore and the input string.
    """
    return str(datetime.date.today().isoformat()) + "_" + s


def log_message(message):
    """
    Returns a log message string with a timestamp.

    Parameters:
        message (str): The message to log.

    Returns:
        str: A string combining the current timestamp and the message.
    """
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    # Return the message with the timestamp
    return f"{timestamp} - {message}"


def check_output_directory(path_to_folder):
    """
    Checks if a directory exists and creates it if it does not.

    Parameters:
        path_to_folder (str or Path): The path to the folder.

    Returns:
        None
    """
    if not os.path.exists(path_to_folder):
        os.makedirs(path_to_folder)
    return


def sort_dict_by_value_len(dict_):
    """
    Returns a list of dictionary items sorted by the length of their values in descending order.

    Parameters:
        dict_ (dict): The dictionary to sort.

    Returns:
        list: A list of tuples (key, value) sorted by the length of the value, descending.
    """
    return sorted(dict_.items(), key=lambda kv: (len(kv[1]), kv[0]))[::-1]


def reverse_complement(seq):
    """
    Returns the reverse complement of a DNA sequence.

    Parameters:
        seq (str): DNA sequence.

    Returns:
        str: Reverse complement of the sequence.
    """
    translate = {'A': 'T ', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N',
                 'S': 'S', 'W': 'W', 'M': 'K', 'K': 'M', 'R': 'Y',
                 'Y': 'R', 'B': 'V', 'V': 'B', 'D': 'H', 'H': 'D'}
    return ''.join([translate[n] for n in seq[::-1]])


def translate_well_to_2d_array(well, coding):
    """
    Translates a well identifier into 2D coordinates for a 96-well plate.

    Parameters:
        well (str): The well identifier (e.g., "C12").
        coding (str): The coding scheme to use ("CK" or "PFS/SE").

    Returns:
        tuple: A tuple (row, column) with zero-indexed row and column numbers.
               For example, "C12" returns (2, 11).
    """
    if (coding == "CK") or (coding == "UNI_v2"):
        translate = {'A': 0, 'B': 1, 'C': 2, 'D': 3, 'E': 4, 'F': 5, 'G': 6, 'H': 7}
        return translate[well[0]], int(well.split(well[0])[1]) - 1


def translate_nt_to_aa(seq):
    """
    Translates a DNA sequence into an amino acid sequence using the first reading frame.

    Parameters:
        seq (str): DNA sequence.

    Returns:
        str: The amino acid sequence using single-letter codes.
             Returns "N/A" if the input sequence is "N/A".
    """
    codon_table = {
        'ATA': 'I', 'ATC': 'I', 'ATT': 'I', 'ATG': 'M',
        'ACA': 'T', 'ACC': 'T', 'ACG': 'T', 'ACT': 'T',
        'AAC': 'N', 'AAT': 'N', 'AAA': 'K', 'AAG': 'K',
        'AGC': 'S', 'AGT': 'S', 'AGA': 'R', 'AGG': 'R',
        'CTA': 'L', 'CTC': 'L', 'CTG': 'L', 'CTT': 'L',
        'CCA': 'P', 'CCC': 'P', 'CCG': 'P', 'CCT': 'P',
        'CAC': 'H', 'CAT': 'H', 'CAA': 'Q', 'CAG': 'Q',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R',
        'GTA': 'V', 'GTC': 'V', 'GTG': 'V', 'GTT': 'V',
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A',
        'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'GGA': 'G', 'GGC': 'G', 'GGG': 'G', 'GGT': 'G',
        'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TTC': 'F', 'TTT': 'F', 'TTA': 'L', 'TTG': 'L',
        'TAC': 'Y', 'TAT': 'Y', 'TAA': '*', 'TAG': '*',
        'TGC': 'C', 'TGT': 'C', 'TGA': '*', 'TGG': 'W',
    }
    dna = str(seq)
    protein_sequence = ""
    if dna == "N/A":
        protein_sequence = "N/A"
    else:
        for codon_start in range(0, len(dna), 3):
            protein_sequence += codon_table.get(dna[codon_start:codon_start + 3].upper(), "?")
    return protein_sequence.rstrip('?')


def hamming(s1, s2):
    """
    Calculates the Hamming distance between two strings.

    Parameters:
        s1 (str): The first string.
        s2 (str): The second string.

    Returns:
        int: The Hamming distance if the strings are of equal length.
             If the strings are of different lengths, prints an error message and exits.
    """
    if len(s1) == len(s2):
        hadi = 0
        for idx, letter in enumerate(s1):
            if letter != s2[idx]:
                hadi += 1
        return hadi
    else:
        print("Hamming distance not computable due to differences in sequence lengths.")
        return sys.exit(2)


def chunk(l, n):
    """
    Yields successive n-sized chunks from a list.

    Parameters:
        l (list): The list to chunk.
        n (int): The size of each chunk.

    Yields:
        list: A chunk of the list with n elements.
    """
    for i in range(0, len(l), n):
        yield l[i:i + n]


def get_sequence(start, end, sequence, frame_shift):
    """
    Returns a substring of a sequence, applying a frame shift.

    The substring is extracted from (start + frame_shift) to (end + frame_shift + 1).

    Parameters:
        start (int): The starting index.
        end (int): The ending index.
        sequence (str): The sequence from which to extract.
        frame_shift (int): The frame shift to apply.

    Returns:
        str: The extracted substring, or "N/A" if start or end is a string.
    """
    if any([type(x) == str for x in [start, end]]):
        return "N/A"
    else:
        start, end, frame_shift = int(start), int(end), int(frame_shift)
        return sequence[start + frame_shift:end + frame_shift + 1]


def decompress(btop):
    """
    Decompresses a BTOP string from IgBLAST.

    In the current version, numeric parts are replaced by repeated ".." indicating identical alignment pairs.
    Returns a list of strings corresponding to each nucleotide position.

    Parameters:
        btop (str): The BTOP string.

    Returns:
        list: A decompressed list representing the BTOP information.
    """
    btop_list = ["".join(x) for _, x in itertools.groupby(btop, key=str.isdigit)]
    decompressed_btop = []
    for part in btop_list:
        if part.isdigit():
            decompressed_btop += [".." for _ in range(int(part))]
        else:
            decompressed_btop.append(part)
    return decompressed_btop


def get_indels(btop, decompression=True):
    """
    Determines the number of insertions and deletions in a BTOP string.

    The function decompresses the BTOP string (if decompression is True) and then counts:
      - Insertions: pairs ending with '-'
      - Deletions: pairs starting with '-'

    Parameters:
        btop (str): The BTOP string.
        decompression (bool): If True, decompresses the BTOP string before processing (default is True).

    Returns:
        tuple: A tuple (insertions, deletions) with the counts of insertions and deletions.
    """
    insertions = 0
    deletions = 0

    if decompression:
        btop = decompress(btop)
    for part in btop:
        mms = list(chunk(part, 2))  # Chunk the BTOP part into pairs.
        for mm in mms:
            if mm.endswith("-"):
                # Insertion detected.
                insertions += 1
            elif mm.startswith("-"):
                # Deletion detected.
                deletions += 1

    return insertions, deletions

def extract_infos_from_file(ab1_file):
    """
    Extracts information from the .ab1 file name (stem).

    The file name is expected to follow a specific format:
    COHORT_SUBJECT_TIME_POINT_TISSUE_SUBSET_PLATE_WELL_CHAIN_PRIMERSET_SOURCE_REVPRIMER_SEQCOUNTER
    For example:
      "HCMV-Control_Mouse-473_t0_Spleen_gB_1_A1_HC_oPR_2ND_IgGInt_EF_SEQ-1.ab1"

    Parameters:
        ab1_file (Path): A Path object representing the .ab1 file.

    Returns:
        tuple: A tuple containing:
               (cohort, subject, time_point, tissue, bait, plate, well, primer_set, chain,
                source, subsource, rev_primer, seq_company, seq_counter)
    """
    file_name = ab1_file.stem
    cohort = file_name.split("_")[0]
    subject = file_name.split("_")[1]
    time_point = file_name.split("_")[2]
    tissue = file_name.split("_")[3]
    bait = file_name.split("_")[4]
    plate = int(file_name.split("_")[5])
    well = file_name.split("_")[6]
    chain = file_name.split("_")[7]
    primer_set = file_name.split("_")[8]
    source = file_name.split("_")[9]
    if "-" in source:
        subsource = "-".join(source.split("-")[1:])
        source = source.split("-")[0]
    else:
        subsource = "N/A"
    rev_primer = file_name.split("_")[10]
    seq_company = file_name.split("_")[11]
    seq_counter = int(file_name.split("_")[12].split('-')[1])

    return (cohort, subject, time_point, tissue, bait, plate, well, primer_set, chain,
            source, subsource, rev_primer, seq_company, seq_counter)


def check_abi_file(ab1_file, handle, q_cut_off, mean_q_cut_off, min_length, inner_n, igblast_results):
    """
    Extracts important data from an ab1 file and its corresponding IgBLAST result, returning a list of quality metrics and sequences.

    The function performs the following steps:
      1. Extracts file name information.
      2. Retrieves relevant IgBLAST result information for quality assessment.
      3. Reads the ab1 file to obtain the raw sequence and quality scores.
      4. Trims the sequence and quality scores based on FWR1 start and J gene end.
      5. Calculates the mean quality of the trimmed sequence.
      6. Masks low-quality bases in the trimmed sequence.
      7. Evaluates if all quality criteria are passed.

    Parameters:
        ab1_file (Path): The .ab1 file path.
        handle: File handle for reading the .ab1 file.
        q_cut_off (int): Quality score cutoff for masking low-quality bases.
        mean_q_cut_off (int): Mean quality score cutoff.
        min_length (int): Minimum required length of the trimmed sequence.
        inner_n (int): Maximum allowed number of masked bases ("N") in the inner region.
        igblast_results (dict): Dictionary containing IgBLAST results keyed by sample name.

    Returns:
        list: A list of extracted values including:
              [cohort, subject, time_point, tissue, bait, plate, well, primer_set, chain, source, subsource, rev_primer,
               seq_company, seq_counter, raw_length, raw_mean_quality, igblast_pass, igblast_v_align, trimmed_length,
               trimmed_length_pass, mean_quality, mean_quality_pass, total_N, inner_n_pass, passed, stop_codon, frame,
               productive, trimmed_seq, masked_sequence, original_sequence, ab1_file.name]
    """
    # Initialize default values.
    trimmed_seq = "N/A"
    trimmed_phred_score = "N/A"
    masked_sequence = "N/A"
    mean_quality = "N/A"
    total_N = "N/A"
    mean_quality_pass = np.nan  # Changed on 2024-11-28.
    trimmed_length = "N/A"
    trimmed_length_pass = np.nan  # Changed on 2024-11-28.
    inner_n_pass = np.nan  # Changed on 2024-11-28.
    passed = False

    # Step 1: Extract file name information.
    cohort, subject, time_point, tissue, bait, plate, well, primer_set, chain, source, subsource, rev_primer, seq_company, seq_counter = extract_infos_from_file(ab1_file)
    name = ab1_file.stem

    # Step 2: Retrieve information from the corresponding IgBLAST result.
    orientation = igblast_results[name]['ORIENTATION']
    fwr1_start = igblast_results[name]['FWR1_START']
    j_end = igblast_results[name]['J_END']
    frame = igblast_results[name]['FRAME']
    stop_codon = igblast_results[name]['STOP_CODON']
    productive = igblast_results[name]['PRODUCTIVE']
    igblast_v_align = igblast_results[name]['FULL_V_ALIGNMENT']
    igblast_pass = igblast_results[name]['IGBLAST_PASSED']

    # Step 3: Read the ab1 file and extract the sequence and quality scores.
    for record in SeqIO.parse(handle, "abi"):
        # Get sequence in the correct orientation and remove all whitespaces.
        seq = reverse_complement(str(record.seq)).replace(" ", "") if orientation == '-' else str(record.seq).replace(" ", "")
        raw_length = len(seq)
        phred_score = record.letter_annotations['phred_quality']
        if orientation == '-':
            phred_score.reverse()
        raw_mean_quality = round(sum(phred_score) / len(seq))

    # Proceed only if IgBLAST results indicate a pass.
    if igblast_pass:
        mean_quality_pass = False
        inner_n_pass = False

        # Step 4: Trim the sequence and quality scores between FWR1 start and J gene end.
        trimmed_seq = seq[fwr1_start - 1:j_end]  # -1 to convert to 0-indexed.
        trimmed_phred_score = phred_score[fwr1_start - 1:j_end]
        trimmed_length = len(trimmed_seq)
        trimmed_length_pass = True if trimmed_length >= min_length else False

        # Step 5: Calculate mean quality of the trimmed sequence.
        mean_quality = round(sum(trimmed_phred_score) / len(trimmed_seq))
        if mean_quality >= mean_q_cut_off:
            mean_quality_pass = True

        # Step 6: Mask bases in the trimmed sequence based on quality cutoff.
        masked_sequence = ""
        for position, base in enumerate(trimmed_seq):
            q = trimmed_phred_score[position]
            if q < q_cut_off:
                masked_sequence += "N"
            else:
                masked_sequence += base
        total_N = masked_sequence.count("N")
        if total_N <= inner_n:
            inner_n_pass = True

        # Step 7: Check if all quality criteria are passed.
        passed = mean_quality_pass and trimmed_length_pass and inner_n_pass

    return [
        cohort, subject, time_point, tissue, bait, plate, well, primer_set, chain, source, subsource, rev_primer,
        seq_company, seq_counter, raw_length, raw_mean_quality, igblast_pass, igblast_v_align, trimmed_length,
        trimmed_length_pass, mean_quality, mean_quality_pass, total_N, inner_n_pass, passed, stop_codon, frame,
        productive, trimmed_seq, masked_sequence, seq, ab1_file.name
    ]


def generate_q_check_df(data_dir, q_co, mean_q_co, min_l, in_n, igblast_d):
    """
    Walks through a directory to process each .ab1 file, performing quality checks, and returns a sorted DataFrame.

    The function iterates over all .ab1 files in the specified directory, applies the check_abi_file function to each,
    collects the resulting data into a list, converts the list into a DataFrame with predefined columns, and sorts the DataFrame
    by SUBJECT, TIME_POINT, SUBSET, PLATE, and WELL.

    Parameters:
        data_dir (str or Path): Directory containing .ab1 files.
        q_co (int): Quality cutoff value.
        mean_q_co (int): Mean quality cutoff.
        min_l (int): Minimum required trimmed sequence length.
        in_n (int): Maximum allowed number of masked bases ("N") in the inner region.
        igblast_d (dict): Dictionary of IgBLAST results keyed by sample name.

    Returns:
        pd.DataFrame: A sorted DataFrame containing the quality check information.
    """
    q_check_list = []
    ab1_files = list(Path(data_dir).rglob("*.ab1"))
    # Process each .ab1 file in the directory.
    for ab1_file in ab1_files:
        with ab1_file.open("rb") as input_file:
            q_check_list.append(
                check_abi_file(ab1_file, input_file, q_co, mean_q_co, min_l, in_n, igblast_d)
            )
    columns = [
        "COHORT", "SUBJECT", "TIME_POINT", "TISSUE", "SUBSET", "PLATE", "WELL", "PRIMER_SET", "CHAIN_PCR",
        "SOURCE", "SUBSOURCE", "REV_PRIMER", "COMPANY", "SEQ_COUNTER", "RAW_LENGTH", "RAW_MEAN_Q",
        "IGBLAST_PASSED", "FULL_V_ALIGNMENT", "V_to_J_LENGTH", "LENGTH_PASSED", "MEAN_Q", "MEAN_PHRED_PASSED", "INNER_N",
        "INNER_N_PASSED", "QCHECK_PASSED", "STOP_CODON", "FRAME", "PRODUCTIVE", "TRIMMED_SEQ", "MASKED_SEQ",
        "ORIGINAL_SEQ", "FILE_NAME"
    ]
    return pd.DataFrame(q_check_list, columns=columns).sort_values(
        by=['SUBJECT', 'TIME_POINT', 'SUBSET', 'PLATE', 'WELL']
    )


def highlight_passed(row):
    """
    Helper function for formatting Excel rows based on quality check results.

    If QCHECK_PASSED is True and PRODUCTIVE is "Yes", the row is colored green.
    If QCHECK_PASSED is True and PRODUCTIVE is "No", the row is colored orange.
    Otherwise, the row is colored red.

    Parameters:
        row (pd.Series): A row from the DataFrame.

    Returns:
        list: A list of style strings for each cell in the row.
    """
    if row['QCHECK_PASSED'] and (row['PRODUCTIVE'] == "Yes"):
        return ['background-color: ' + green for _ in row]
    elif row['QCHECK_PASSED'] and (row['PRODUCTIVE'] == "No"):
        return ['background-color: ' + orange for _ in row]
    else:
        return ['background-color: ' + red for _ in row]


def make_quality_check_excel(q_check_df, path_to_excel_file):
    """
    Generates an Excel file from the quality check DataFrame with conditional formatting.

    The function applies color-coding to rows based on quality check results and exports the styled DataFrame
    to an Excel file at the specified path.

    Parameters:
        q_check_df (pd.DataFrame): DataFrame containing quality check information.
        path_to_excel_file (str or Path): The file path where the Excel file will be saved.

    Returns:
        tuple: (True, message) if export is successful; (False, error message) if an error occurs.
    """
    try:
        output_style = q_check_df.reset_index(drop=True).style.apply(highlight_passed, axis=1)
        with pd.ExcelWriter(path_to_excel_file) as writer:
            output_style.to_excel(writer, sheet_name='Summary')
        return True, "Quality check Excel file exported"
    except Exception as e:
        return False, f"Error in make_quality_check_excel: {e}"

# ===================
# PLOTTING FUNCTIONS
# ===================

def plot_pie_chart(
        pie_sizes,
        pie_colors,
        center_text="",
        center_fontsize=20,
        counterclockwise=False,
        edges=None,
        legend=None,
        bbox_to_anchor_tuple=(0.5, 0),
        figsize=(4, 4),
        wedge_width=0.5,
        legend_font_size=12,
        max_legend_items_per_col=8,
        fixed_pie_position=False,
        adjust_kwargs=None
):
    """
    Plots a donut (or simple pie) chart with custom parameters and returns a Matplotlib Figure.

    This function creates a pie chart (donut chart if wedge_width is set) with the specified slice sizes and colors.
    It allows for center text, custom edge properties, and an optional legend. The figure size and layout can be
    adjusted via parameters.

    Parameters:
        pie_sizes (list): Numeric values for each slice of the pie.
        pie_colors (list): Color codes for each slice (e.g., ['#FF0000', '#00FF00', ...]).
        center_text (str): Text to display in the center of the donut (default is an empty string).
        center_fontsize (int): Font size for the center text (default is 20).
        counterclockwise (bool): If True, the slices are drawn counterclockwise (Matplotlib uses counterclockwise by default).
        edges (dict): Dictionary of edge properties, e.g. {'color': 'black', 'ewidth': 0.5}.
        legend (dict): Dictionary with legend parameters, e.g. {'labels': [...], 'loc': 'upper left', 'title': 'Legend'}.
        bbox_to_anchor_tuple (tuple): Tuple for legend offset (default is (0.5, 0)).
        figsize (tuple): Figure size in inches, e.g. (4, 4) (default is (4, 4)).
        wedge_width (float): Width of the donut ring (default is 0.5; if set to 1, produces a full pie).
        legend_font_size (int): Font size for the legend text (default is 12).
        max_legend_items_per_col (int): Maximum number of legend items per column (default is 8).
        fixed_pie_position (bool or dict): If provided, manually sets the pie chart's axis position.
        adjust_kwargs (dict): Optional dictionary with parameters for fig.subplots_adjust (default is {'bottom': 0.3}).

    Returns:
        matplotlib.figure.Figure: The generated Matplotlib Figure object.
    """
    # Create figure and axis.
    fig, ax = plt.subplots(figsize=figsize)

    # Set default wedge properties if edges are provided.
    wedgeprops = {}
    if edges is not None:
        wedgeprops['edgecolor'] = edges.get('color', 'black')
        wedgeprops['linewidth'] = edges.get('ewidth', 0.5)
        # Extend wedge properties here if needed.

    # Plot the pie chart with the specified parameters.
    wedges, texts = ax.pie(
        pie_sizes,
        colors=pie_colors,
        radius=1,
        startangle=90,
        counterclock=counterclockwise,
        wedgeprops=dict(width=wedge_width, **wedgeprops)
    )

    # Ensure the pie chart is drawn as a donut by setting the aspect ratio to equal.
    ax.axis('equal')

    # Add center text if provided.
    if center_text:
        ax.text(0.5, 0.5, str(center_text),
                horizontalalignment='center',
                verticalalignment='center',
                fontsize=center_fontsize, color='black',
                transform=ax.transAxes)

    # Adjust the figure layout using default or provided adjustments.
    base_adjust = {"bottom": 0.3}
    if adjust_kwargs is not None:
        base_adjust.update(adjust_kwargs)
    fig.subplots_adjust(**base_adjust)

    # Manually adjust the axis position if a fixed pie position is provided.
    if fixed_pie_position:
        ax.set_position(fixed_pie_position)

    # Add legend if specified.
    if legend is not None:
        labels = legend.get('labels', None)
        # Default legend placement parameters.
        loc = legend.get('loc', 'upper center')
        bbox_to_anchor = legend.get('bbox_to_anchor', bbox_to_anchor_tuple)
        title = legend.get('title', None)

        if labels:
            n = len(labels)
            # Determine the number of columns: 1 column for up to max_legend_items_per_col labels, etc.
            ncol = math.ceil(n / max_legend_items_per_col)
            ax.legend(wedges, labels, loc=loc, bbox_to_anchor=bbox_to_anchor,
                      title=title, fontsize=legend_font_size, ncol=ncol)

    return fig


def plot_qc_statistics(q_check_df):
    """
    Identifies the reasons for QC failure and returns a horizontal bar plot.

    This function filters the quality check DataFrame for rows that did not pass QC, then counts
    the number of failures for each of the following criteria:
      - IgBLAST failure
      - Alignment too short (LENGTH_PASSED)
      - Poor overall quality (MEAN_PHRED_PASSED)
      - Excessive uncertain bases (INNER_N_PASSED)

    A horizontal bar chart is then generated displaying the count of sequences that failed each criterion.

    Parameters:
        q_check_df (pd.DataFrame): DataFrame containing quality check information, including columns
                                   "QCHECK_PASSED", "MEAN_PHRED_PASSED", "LENGTH_PASSED", "INNER_N_PASSED",
                                   and "IGBLAST_PASSED".

    Returns:
        matplotlib.figure.Figure: A Matplotlib Figure object containing the bar plot.
    """
    # Filter the DataFrame for sequences that did not pass quality check.
    np_df = q_check_df[~q_check_df["QCHECK_PASSED"]]

    # Count failures for each QC criterion.
    qp = len(np_df[~np_df["MEAN_PHRED_PASSED"].astype(bool)])
    ilp = len(np_df[~np_df["LENGTH_PASSED"].astype(bool)])
    inp = len(np_df[~np_df["INNER_N_PASSED"].astype(bool)])
    igbp = len(np_df[~np_df["IGBLAST_PASSED"].astype(bool)])

    # List of failure counts and corresponding labels.
    fail_reason_list = [igbp, ilp, qp, inp]
    fail_names = ["IgBLAST\nfailed", "Alignment\ntoo short", "Poor overall\nquality", "Excessive\nuncertain bases"]

    # Create a horizontal bar plot.
    fig, ax = plt.subplots(figsize=(4, 4))
    ax.barh(range(len(fail_reason_list)),
            fail_reason_list,
            color=[light_grey, red, red, red],
            edgecolor='black',
            linewidth=0.5)
    ax.yaxis.set_ticks(range(len(fail_names)))
    ax.set_yticklabels(fail_names)
    ax.set_xlabel("Number of sequences")

    return fig

# TODO: should be moved to analyze_ab1_files.py
def plot_qc_donut_chart(df):
    """
    Computes QC statistics from the DataFrame and generates a donut chart summarizing the results.

    This function calculates counts for the following categories:
      - "QC passed and productive": Sequences that passed quality check and are productive.
      - "QC passed but non-productive": Sequences that passed quality check but are non-productive.
      - "QC not passed and non-productive": Sequences that did not pass quality check and are non-productive.
      - "IgBLAST failed": Sequences that failed the IgBLAST step (computed as the total number minus the sum of the other counts).

    The counts are used to create a donut chart via the plot_pie_chart function. A legend is generated showing the category, count,
    and percentage of total sequences.

    Parameters:
        df (pd.DataFrame): DataFrame containing QC results. Expected to include the columns 'QCHECK_PASSED' and 'PRODUCTIVE'.

    Returns:
        matplotlib.figure.Figure: The generated donut chart as a Matplotlib Figure object.
    """
    total = len(df)

    # Calculate counts for each QC category.
    category_counts = {
        "QC passed and productive": len(df[(df['QCHECK_PASSED'] == True) & (df['PRODUCTIVE'] == "Yes")]),
        "QC passed but non-productive": len(df[(df['QCHECK_PASSED'] == True) & (df['PRODUCTIVE'] == "No")]),
        "QC not passed and non-productive": len(df[(df['QCHECK_PASSED'] == False) & (df['PRODUCTIVE'] == "No")])
    }
    category_counts["IgBLAST failed"] = total - sum(category_counts.values())

    sizes = list(category_counts.values())
    colors = [green, orange, red, light_grey]

    # Create legend labels with counts and percentages.
    legend_labels = []
    for cat, count in category_counts.items():
        percentage = (count / total * 100) if total > 0 else 0
        legend_labels.append(f"{cat}: {count} ({round(percentage, 1)}%)")

    legend_dict = {
        'labels': legend_labels,
        'loc': 'upper center',
        'title': 'Legend'
    }

    # Generate the donut chart using the plot_pie_chart function.
    fig = plot_pie_chart(
        pie_sizes=sizes,
        pie_colors=colors,
        center_text=str(total),
        counterclockwise=False,
        edges={'color': 'black', 'ewidth': 0.5},
        legend=legend_dict,
        figsize=(4, 4),
        wedge_width=0.5
    )

    return fig

### BLAST FUNCTIONS

def perform_blast(query, db, strand="plus", evalue="20", word_size="9", penalty="-1", gapopen="4",
                  gapextend="1", dust="no", outfmt="'7 sseqid mismatch btop evalue score'", max_target_seqs="3",
                  culling_limit="3"):
    """
    Executes a BLAST search using blastn with specified parameters and returns the parsed output.

    This function constructs a command-line BLAST query with the provided parameters (such as database,
    E-value, word size, penalty, gap penalties, dust filtering, etc.) and executes it via a subprocess.
    The BLAST output is read from the standard output and split by the marker "# BLASTN", returning the
    result entries following this marker.

    Parameters:
        query (str or Path): Path to the FASTA file containing the query sequences.
        db (str): The BLAST database to search against.
        strand (str): Strand to search ("plus" by default).
        evalue (str): E-value threshold (default: "20").
        word_size (str): Word size parameter (default: "9").
        penalty (str): Penalty for mismatches (default: "-1").
        gapopen (str): Gap open penalty (default: "4").
        gapextend (str): Gap extension penalty (default: "1").
        dust (str): Dust filtering option (default: "no").
        outfmt (str): Output format for BLAST results (default: "'7 sseqid mismatch btop evalue score'").
        max_target_seqs (str): Maximum number of target sequences to return (default: "3").
        culling_limit (str): Culling limit for BLAST results (default: "3").

    Returns:
        list: A list of BLAST result entries extracted by splitting the output using "# BLASTN" as the delimiter.
    """
    blast_command = " ".join(
        ["blastn",
         "-db", db,
         "-strand", strand,
         "-evalue", evalue,
         "-word_size", word_size,
         "-penalty", penalty,
         "-gapopen", gapopen,
         "-gapextend", gapextend,
         "-dust", dust,
         "-outfmt", outfmt,
         "-max_target_seqs", max_target_seqs,
         "-culling_limit", culling_limit,
         "-query", str(query)])

    print(log_message("Running BLAST"), flush=True)

    blast = Popen(blast_command, stdout=PIPE, stderr=PIPE, shell=True)

    bl_results = blast.stdout.read().decode()
    bl_error = blast.stderr.read().decode()
    # Optionally print error output: print("stderr from blast:\n", bl_error, flush=True)
    return bl_results.split('# BLASTN')[1:]


def get_isotype_with_blast(input_column, output_folder, blast_db="CH1_DB"):
    """
    Determines isotype information by performing a BLAST search on sequences from a Pandas Series.

    This function converts an input Pandas Series (where each value represents a sequence) into a FASTA file.
    The FASTA file is saved to the specified output folder and then blasted against the provided BLAST database
    using the perform_blast function. The function parses the BLAST output and returns a list containing the top hit(s)
    (isotype annotations) for each query sequence.

    Parameters:
        input_column (pd.Series): Pandas Series containing sequences to be queried. Each sequence is labeled with its index.
        output_folder (str or Path): The folder in which to save the temporary FASTA file.
        blast_db (str): The BLAST database to search against (default: "CH1_DB").

    Returns:
        list: A list of isotype annotations corresponding to each sequence in the input column.
    """
    # Export the series to a temporary FASTA file.
    output = []
    for k, v in input_column.items():
        output.append(">" + str(k) + "\n" + str(v))
    path_to_blast_file = Path(output_folder).joinpath(date_stamp("blast-query.fasta"))
    with open(path_to_blast_file, "w") as handle:
        handle.write("\n".join(output))

    # Execute BLAST search.
    blast_output = perform_blast(path_to_blast_file, blast_db)
    isotype_list = []
    for entry in blast_output:
        if entry.split("\n")[3].startswith("# Fields:"):
            number_of_hits = int(entry.split("\n")[4].split(" ")[1])
            if number_of_hits == 0:
                isotype_list.append("N/A")
            else:
                min_e = 20
                all_hits = []
                for i in range(number_of_hits):
                    e = float(entry.split("\n")[5 + i].split("\t")[3])
                    hit = entry.split("\n")[5 + i].split("\t")[0]
                    if e <= min_e:
                        all_hits.append(hit)
                        min_e = e
                isotype_list.append(", ".join(all_hits))
        else:
            isotype_list.append("N/A")
    return isotype_list


def dict_to_str(d, indent=0):
    """
    Recursively converts a dictionary into a list of strings with indentation.

    Each key-value pair is converted into a string. If a value is a dictionary, the function calls itself
    recursively, increasing the indentation level. The output is a list of strings representing the structure
    of the dictionary.

    Parameters:
        d (dict): The dictionary to convert.
        indent (int): The current indentation level (number of spaces). Default is 0.

    Returns:
        list: A list of strings representing the dictionary with appropriate indentation.
    """
    lines = []
    for key, value in d.items():
        prefix = " " * indent + f"{key}: "
        if isinstance(value, dict):
            lines.append(prefix)
            lines.extend(dict_to_str(value, indent=indent + 2))
        else:
            lines.append(prefix + f"{value}")
    return lines


def write_clustering_log(log_file_path, settings=False, input_files=False, output_files=False):
    """
    Creates a clustering log file containing a timestamp, active settings, input file names, and exported file names.

    This function writes a log file at the specified path, which includes:
      - A timestamp in the format YYMMDD-HHMMSS.
      - A list of input file names.
      - A recursively formatted listing of active settings.
      - A list of exported file names.

    Parameters:
        log_file_path (str or Path): The path to the log file to be created.
        settings (dict, optional): A dictionary with active settings (can be nested). Defaults to False.
        input_files (list, optional): A list of input file names (strings). Defaults to False.
        output_files (list, optional): A list of output file names (strings). Defaults to False.

    Returns:
        None
    """
    # Create a timestamp string in the format YYMMDD-HHMMSS.
    timestamp = datetime.datetime.now().strftime("%y%m%d-%H%M%S")

    with open(log_file_path, "w") as log_file:
        log_file.write(f"TIMESTAMP: {timestamp}\n\n")

        # Write input file names.
        log_file.write("Input files:\n")
        if input_files:
            for infile in input_files:
                log_file.write(f"{infile}\n")
        else:
            log_file.write("- No input files passed -\n")
        log_file.write("\n")

        # Write active settings.
        log_file.write("Active filter:\n")
        if settings:
            settings_lines = dict_to_str(settings, indent=0)
            log_file.write("\n".join(settings_lines))
        else:
            log_file.write("- No settings passed -")
        log_file.write("\n\n")

        # Write exported file names.
        log_file.write("Exported files:\n")
        if output_files:
            for outfile in output_files:
                log_file.write(f"{outfile}\n")
        else:
            log_file.write("- No output files passed -\n")