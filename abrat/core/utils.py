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
    """Wendet die globalen Theme-Konfigurationen für Matplotlib, Altair und Plotly an."""
    # Matplotlib
    plt.rcParams['font.family'] = 'DejaVu Sans'
    plt.rcParams['font.size'] = 12
    plt.rcParams['lines.linewidth'] = 0.5
    plt.rcParams['axes.linewidth'] = 0.5
    plt.rcParams['xtick.major.width'] = 0.5
    plt.rcParams['ytick.major.width'] = 0.5
    plt.rcParams['xtick.minor.width'] = 0.5
    plt.rcParams['ytick.minor.width'] = 0.5

    # Altair: definiere und aktiviere das Theme
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

    # Plotly: definiere und setze ein Template als Standard
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

    :param original_filename: Original filename (e.g., "myfile.xlsx")
    :return: New filename with the timestamp inserted (e.g., "myfile_clustered-230308-153045.xlsx")
    """
    base, ext = os.path.splitext(original_filename)
    timestamp = datetime.datetime.now().strftime("-%y%m%d-%H%M%S")
    return f"{base}_{insert_text}{timestamp}{ext}"

def date_stamp(s):
    return str(datetime.date.today().isoformat())+"_"+s

def log_message(message):
    # for printing a message with a timestemp
    timestamp = datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    # Nachricht mit Zeitstempel ausgeben
    return(f"{timestamp} - {message}")

def check_output_directory(path_to_folder):
    """Checks if dir is present and generates it if not."""
    if not os.path.exists(path_to_folder):
        os.makedirs(path_to_folder)
    return

def sort_dict_by_value_len(dict_):
    """Takes dict_ and returns sorted dict by lengths of values."""
    return sorted(dict_.items(), key=lambda kv: (len(kv[1]), kv[0]))[::-1]

def reverse_complement(seq):
    """Returns DNA reverse complement of seq."""
    translate = {'A': 'T ', 'T': 'A', 'G': 'C', 'C': 'G', 'N': 'N',
                 'S': 'S', 'W': 'W', 'M': 'K', 'K': 'M', 'R': 'Y',
                 'Y': 'R', 'B': 'V', 'V': 'B', 'D': 'H', 'H': 'D'}
    return ''.join([translate[n] for n in seq[::-1]])

def translate_well_to_2d_array(well, coding):
    """Takes a coding (either CK or PFS/SE) and translates the well numbering
    into 2D coordinates from a 96 well plate. And returns tuple with
    row and column. E.g. 36/C12 => (3, 12)
    """

    if (coding == "CK") or (coding == "UNI_v2"):
        translate = {'A': 0, 'B': 1, 'C': 2, 'D': 3, 'E': 4, 'F': 5, 'G': 6, 'H': 7}
        return translate[well[0]], int(well.split(well[0])[1]) - 1

def translate_nt_to_aa(seq):
    """Takes DNA sequence in seq and returns the first frame aminoacid
    translation.
    """

    # Translation table
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
            protein_sequence += codon_table.get(dna[codon_start:codon_start + 3].upper(), "?")  # upper case changed
    return protein_sequence.rstrip('?')

def hamming(s1, s2):
    """Calculates Hamming distance"""

    if len(s1) == len(s2):
        hadi = 0
        for idx, letter in enumerate(s1):
            if letter != s2[idx]:
                hadi += 1
        return hadi
    else:
        print("Hamming distance not computable, due to differences in sequence lengths.")
        return sys.exit(2)

def chunk(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]

def get_sequence(start, end, sequence, frame_shift):
    """ Returns sequence from start to end, including end"""
    if any([type(x) == str for x in [start, end]]):
        return "N/A"
    else:
        start, end, frame_shift = int(start), int(end), int(frame_shift)
        return sequence[start + frame_shift:end + frame_shift + 1]

# BTOP Functions #
def decompress(btop):
    """Function to decompress BTOP string from igblast.
    In current version simply adds ".." for identical alignment pair.
    Returns a list of pairs which correspond to each nt position.
    """
    btop_list = ["".join(x) for _, x in itertools.groupby(btop, key=str.isdigit)]
    # go through BTOP list and add mismatch_count to reference_dict
    decompressed_btop = []
    for part in btop_list:
        if part.isdigit():
            decompressed_btop += [".." for _ in range(int(part))]
        else:
            decompressed_btop.append(part)
    return decompressed_btop

def get_indels(btop, decompression=True):
    """Function reads a BTOP and decompresses it.
    Set decompression = False if you want to pass a decompressed BTOP directly.
    will return number of insertions, deletions.
    """
    insertions = 0
    deletions = 0

    if decompression:
        btop = decompress(btop)
    for part in btop:
        mms = list(chunk(part, 2))  # chunk btop in packages of 2 and decipher according to character
        for mm in mms:
            if mm.endswith("-"):
                # => insertion
                insertions += 1
            elif mm.startswith("-"):
                # => deletion
                deletions += 1

    return insertions, deletions

def extract_infos_from_file(ab1_file):
    """ Function to return infos from ab1_file file_name (stem)."""
    # COHORT_PATIENT_ZEITPUNKT_TISSUE_SUBSET_PLATTE_WELL_KETTE_PRIMERSET_NUKLEINSÄUREQUELLE_SEQPRIMER_COMPANY_SEQ-x
    # HCMV-Control_Mouse-473_t0_Spleen_gB_1_A1_HC_oPR_2ND_IgGInt_EF_SEQ-1.ab1
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

    return cohort, subject, time_point, tissue, bait, plate, well, primer_set, chain, source, subsource, \
        rev_primer, seq_company, seq_counter

def check_abi_file(ab1_file, handle, q_cut_off, mean_q_cut_off, min_length, inner_n, igblast_results):
    """Function to extract important data from igblast/abi file and return masked sequence and infos."""

    # initialize values
    trimmed_seq = "N/A"
    trimmed_phred_score = "N/A"
    masked_sequence = "N/A"
    mean_quality = "N/A"
    total_N = "N/A"
    mean_quality_pass = np.nan # edit 2024-11-28
    trimmed_length = "N/A"
    trimmed_length_pass = np.nan # edit 2024-11-28, was "N/A" before
    inner_n_pass = np.nan # edit 2024-11-28, was "N/A" before
    passed = False

    # Step 1: Get infos from ab1-file name
    cohort, subject, time_point, tissue, bait, plate, well, primer_set, chain, source, subsource, rev_primer, \
        seq_company, seq_counter = extract_infos_from_file(ab1_file)
    name = ab1_file.stem

    # Step 2: Get infos from corresponding Blast result and make first quality check
    orientation = igblast_results[name]['ORIENTATION']
    fwr1_start = igblast_results[name]['FWR1_START']
    j_end = igblast_results[name]['J_END']
    frame = igblast_results[name]['FRAME']
    stop_codon = igblast_results[name]['STOP_CODON']
    productive = igblast_results[name]['PRODUCTIVE']
    igblast_v_align = igblast_results[name]['FULL_V_ALIGNMENT']
    igblast_pass = igblast_results[name]['IGBLAST_PASSED']

    # Step 3: Get data from ab1 file
    for record in SeqIO.parse(handle, "abi"):
        # get sequence in correct orientation and remove all whitespaces (=> latter was not neccesary before?!)
        seq = reverse_complement(str(record.seq)).replace(" ", "") if orientation == '-' else \
            str(record.seq).replace(" ", "")

        raw_length = len(seq)
        phred_score = record.letter_annotations['phred_quality']
        if orientation == '-':
            phred_score.reverse()
        raw_mean_quality = round(sum(phred_score) / len(seq))
        # n = record.name

    # Continue only, if igblast results make sense #
    if igblast_pass:
        mean_quality_pass = False
        inner_n_pass = False

        # Step 4: trimm seq and phred score to fwr1_start and j_end
        trimmed_seq = seq[fwr1_start - 1:j_end]  # -1 to account for position "0"
        # print(seq, trimmed_seq)
        trimmed_phred_score = phred_score[fwr1_start - 1:j_end]  # could change that to constant region?!
        trimmed_length = len(trimmed_seq)
        trimmed_length_pass = True if trimmed_length >= min_length else False

        # Step 5: Mean quality of trimmed sequence
        mean_quality = round(sum(trimmed_phred_score) / len(trimmed_seq))
        if mean_quality >= mean_q_cut_off:
            mean_quality_pass = True

        # Step 6: Masking
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

        # Step 7: check if all criteria are passed
        passed = mean_quality_pass and trimmed_length_pass and inner_n_pass

    return [cohort, subject, time_point, tissue, bait, plate, well, primer_set, chain, source, subsource, rev_primer,
            seq_company, seq_counter, raw_length, raw_mean_quality, igblast_pass, igblast_v_align, trimmed_length,
            trimmed_length_pass, mean_quality, mean_quality_pass, total_N, inner_n_pass, passed, stop_codon, frame,
            productive, trimmed_seq, masked_sequence, seq, ab1_file.name]


def generate_q_check_df(data_dir, q_co, mean_q_co, min_l, in_n, igblast_d):
    """Walks through data_dir, checks each ab1 file and blast result with check_abi_file, appends the data
    to q_check_list and converts the list of lists to a sorted dataframe with columns.
    """

    q_check_list = []

    ab1_files = list(Path(data_dir).rglob("*.ab1"))
    # go through folder and evaluate all files
    for ab1_file in ab1_files:
        # open file
        with ab1_file.open("rb") as input_file:
            # add values to q_check_list
            q_check_list.append(
                check_abi_file(ab1_file, input_file, q_co, mean_q_co, min_l,
                               in_n, igblast_d))

    # cohort, subject, time_point, bait, plate, well, primer_set, chain, source, rev_primer, seq_company, seq_counter
    columns = ["COHORT", "SUBJECT", "TIME_POINT", "TISSUE", "SUBSET", "PLATE", "WELL", "PRIMER_SET", "CHAIN_PCR",
               "SOURCE", "SUBSOURCE", "REV_PRIMER", "COMPANY", "SEQ_COUNTER", "RAW_LENGTH", "RAW_MEAN_Q",
               "IGBLAST_PASSED", "FULL_V_ALIGNMENT", "V_to_J_LENGTH", "LENGTH_PASSED", "MEAN_Q", "MEAN_PHRED_PASSED", "INNER_N",
               "INNER_N_PASSED", "PASSED", "STOP_CODON", "FRAME", "PRODUCTIVE", "TRIMMED_SEQ", "MASKED_SEQ",
               "ORIGINAL_SEQ", "FILE_NAME"]

    return pd.DataFrame(q_check_list, columns=columns).sort_values(
        by=['SUBJECT', 'TIME_POINT', 'SUBSET', 'PLATE', 'WELL'])

def highlight_passed(row):
    """Helper Function to format excel rows."""
    if row['PASSED'] and (row['PRODUCTIVE'] == "Yes"):
        return ['background-color: ' + green for _ in row]
    elif row['PASSED'] and (row['PRODUCTIVE'] == "No"):
        return ['background-color: ' + orange for _ in row]
    else:
        return ['background-color: ' + red for _ in row]

def make_quality_check_excel(q_check_df, path_to_excel_file):
    """Function to generate an Excel from q_check_df which is
    sorted and colored."""
    try:
        # color if not passed
        output_style = q_check_df.reset_index(drop=True).style.apply(highlight_passed, axis=1)
        with pd.ExcelWriter(path_to_excel_file) as writer:
            output_style.to_excel(writer, sheet_name='Summary')

        return True, f"Quality check excel file exported"
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
        legend_font_size=12 ,
        max_legend_items_per_col=8,
        fixed_pie_position=False,
        adjust_kwargs=None
):
    """
    Plots a donut (or simple pie) chart with the given parameters and returns the figure.

    :param pie_sizes: List of numeric values for each slice of the pie.
    :param pie_colors: List of color codes (e.g. ['#FF0000', '#00FF00', ...]).
    :param center_text: String to display in the center of the donut.
    :param clockwise: If True, the pie is drawn counterclockwise = False. (Matplotlib uses counterclock by default.)
    :param edges: Dictionary with edge props, e.g. {'color':'black', 'ewidth':0.5}.
    :param legend: Dictionary with legend parameters, e.g. {'labels': [...], 'loc': 'upper left'}.
    :param figsize: Tuple specifying figure size, e.g. (4,4).
    :param wedge_width: Width of the donut hole, default=0.5.
    :param adjust_kwargs: Optional dictionary with parameters for fig.subplots_adjust.
                          By default, {'bottom': 0.3} is used.
    :param bbox_to_anchor_tuple=(0.5, 0) to set legend offset
    :return: A Matplotlib Figure object.
    """
    # Create figure and axis
    fig, ax = plt.subplots(figsize=figsize)

    # Default wedge properties
    wedgeprops = {}
    if edges is not None:
        wedgeprops['edgecolor'] = edges.get('color', 'black')
        wedgeprops['linewidth'] = edges.get('ewidth', 0.5)
        # If you want more edge properties, extend here

    # Plot the pie
    wedges, texts = ax.pie(
        pie_sizes,
        colors=pie_colors,
        radius=1,
        startangle=90,
        counterclock=counterclockwise,
        wedgeprops=dict(width=wedge_width, **wedgeprops)
    )

    # Make it look like a donut
    ax.axis('equal')

    # Center text
    if center_text:
        ax.text(0.5, 0.5, str(center_text),
                horizontalalignment='center',
                verticalalignment='center',
                fontsize=center_fontsize, color='black',
                transform=ax.transAxes)

    base_adjust = {"bottom": 0.3} #left=0.2, right=0.8, top=0.8,´
    if adjust_kwargs is not None:
        # Update base_adjust with provided parameters.
        base_adjust.update(adjust_kwargs)

    fig.subplots_adjust(**base_adjust)

    # Manually adjust axis position so that left margin is minimized and the pie chart's area is fixed
    if fixed_pie_position:
        ax.set_position(fixed_pie_position)

    # Legend
    if legend is not None:
        labels = legend.get('labels', None)
        # Set default values for legend placement to be below the chart.
        loc = legend.get('loc', 'upper center')
        bbox_to_anchor = legend.get('bbox_to_anchor', bbox_to_anchor_tuple)
        title = legend.get('title', None)

        if labels:
            n = len(labels)
            # Determine number of columns: 1 column for <=8 labels, 2 columns for 9-16, 3 for 17-24, etc.
            ncol = math.ceil(n / max_legend_items_per_col)
            ax.legend(wedges, labels, loc=loc, bbox_to_anchor=bbox_to_anchor,
                      title=title, fontsize=legend_font_size, ncol=ncol)


    return fig


def plot_qc_statistics(q_check_df):
    """ Identifeis reaseons why QC failed und returns a bar plot."""

    # find reasons for not passing "MEAN_PHRED_PASSED","INIT_LENGTH_PASS", "TRIMMED_LENGTH_PASS", "INNER_N_PASS"
    np_df = q_check_df[~q_check_df["QCHECK_PASSED"]]

    qp = len(np_df[~np_df["MEAN_PHRED_PASSED"].astype(bool)])  # /len(np_df) if len(np_df)>0 else 0 # put 100* in front to get %
    ilp = len(np_df[~np_df["LENGTH_PASSED"].astype(bool)])  # /len(np_df) if len(np_df)>0 else 0
    inp = len(np_df[~np_df["INNER_N_PASSED"].astype(bool)])  # /len(np_df) if len(np_df)>0 else 0
    igbp = len(np_df[~np_df["IGBLAST_PASSED"].astype(bool)])  # /len(np_df) if len(np_df)>0 else 0

    fail_reason_list = [igbp, ilp, qp, inp]
    fail_names = ["IgBLAST\nfailed", "Alignment\nto short", "Bad overall\nquality", "Uncertain\nbases"]

    # Plot summary statistics
    fig, ax = plt.subplots(figsize=(4,4))
    ax.barh([y for y in range(len(fail_reason_list))],
            fail_reason_list,
            color=[light_grey, red, red, red],
            edgecolor='black',
            linewidth=0.5
            )
    ax.yaxis.set_ticks(range(0, len(fail_names), 1))
    ax.set_yticklabels(fail_names)
    #for tick in ax.get_yticklabels():
    #    tick.set_rotation(90)
    ax.set_xlabel("No. of sequences")

    return fig

# TODO: should be moved to analyze_ab1_files.py
def plot_qc_donut_chart(df):
    """
    Computes parameters from the DataFrame and creates a QC donut chart.
    Returns the Matplotlib Figure object. The caller is responsible for saving the figure.

    :param df: DataFrame containing QC results.
    :param output_folder: Path to the folder for saving the figure.
    :return: Matplotlib Figure object.
    """
    total = len(df)

    category_counts = {
        "QC passed and productive": len(df[(df['QCHECK_PASSED'] == True) & (df['PRODUCTIVE'] == "Yes")]),
        "QC passed but non-productive": len(df[(df['QCHECK_PASSED'] == True) & (df['PRODUCTIVE'] == "No")]),
        "QC not passed and non-productive": len(df[(df['QCHECK_PASSED'] == False) & (df['PRODUCTIVE'] == "No")])
    }
    category_counts["IgBLAST failed"] = total - sum(category_counts.values())

    sizes = list(category_counts.values())

    colors = [green, orange, red, light_grey]

    legend_labels = []
    for cat, count in category_counts.items():
        percentage = (count / total * 100) if total > 0 else 0
        legend_labels.append(f"{cat}: {count} ({round(percentage, 1)}%)")

    legend_dict = {
        'labels': legend_labels,
        'loc': 'upper center',
        'title': 'Legend'
    }

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
    """Performs BLAST search on query and returns result."""
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

    # save output to variable
    bl_results = blast.stdout.read().decode()
    bl_error = blast.stderr.read().decode()
    #print("stderr from blast:\n", bl_error, flush=True)
    return bl_results.split('# BLASTN')[1:]

def get_isotype_with_blast(input_column, output_folder, blast_db="CH1_DB"):
    """Functions takes input column (panda series) and converts it to fasta_file.
    fasta file is blasted against blast_db (with perform_blast function)
    and a list with the top hits is returned.
    """
    # export table to temporary file
    output = []
    for k, v in input_column.items():
        output.append(">" + str(k) + "\n" + str(v))
    path_to_blast_file = Path(output_folder).joinpath(date_stamp("blast-query.fasta"))
    with open(path_to_blast_file, "w") as handle:
        handle.write("\n".join(output))

    # blast
    blast_output = perform_blast(path_to_blast_file, blast_db)
    isotype_list = []
    for entry in blast_output:
        #query = entry.split("\n")[1].split(": ")[1]
        if entry.split("\n")[3].startswith("# Fields:"):
            number_of_hits = int(entry.split("\n")[4].split(" ")[1])
            if number_of_hits == 0:
                isotype_list.append("N/A")
            else:
                min_e = 20
                all_hits = []
                for i in range(number_of_hits):
                    e = float(entry.split("\n")[5 + i].split("\t")[3])
                    #score = entry.split("\n")[5 + i].split("\t")[4]
                    hit = entry.split("\n")[5 + i].split("\t")[0]
                    if e <= min_e:
                        all_hits.append(hit)
                        min_e = e
                isotype_list.append(", ".join(all_hits))
        else:
            isotype_list.append("N/A")
    #print(isotype_list)
    return isotype_list


# TODO: change logic for combining dataframes and probably copy to analyze_ab1_filey.py
def combine_ig_blast_and_q_check_data(q_check_df, igblast_dict, combi_ex_name, out_dir, isotype_determination):
    """
    Kombiniert die IgBLAST-Daten mit den Qualitätsprüfungsdaten zu einem DataFrame.
    Optimierte Version: sammelt die einzelnen DataFrames in einer Liste und führt einmalig ein pd.concat aus.
    """
    df_list = []  # Liste für alle Zwischen-DataFrames

    for idx, entry in q_check_df.iterrows():
        # Dateiname ohne Suffix
        file_name = entry['FILE_NAME'].replace('.ab1', '')
        orig_seq = entry['ORIGINAL_SEQ']

        # Hole IgBLAST-Ergebnis (als Dictionary) und extrahiere benötigte Werte
        igblast_result = igblast_dict.get(file_name, {})
        cdr3_nt = igblast_result.get('CDR3_NT', "N/A")
        cdr3_aa = igblast_result.get('CDR3_AA', "N/A")

        # Erstelle ein DataFrame aus dem IgBLAST-Ergebnis
        new_df = pd.DataFrame([igblast_result])
        new_df['SAMPLE_NAME'] = file_name
        new_df['SUBSET'] = entry['SUBSET']
        new_df['TISSUE'] = entry['TISSUE']
        new_df['CHAIN_PCR'] = entry['CHAIN_PCR']
        new_df['IGBLAST_PASSED'] = entry['IGBLAST_PASSED']
        new_df['MEAN_PHRED_PASSED'] = entry['MEAN_PHRED_PASSED']
        new_df['INNER_N_PASSED'] = entry['INNER_N_PASSED']
        new_df['LENGTH_PASSED'] = entry['LENGTH_PASSED']
        new_df['INNER_N'] = entry['INNER_N']
        new_df['QCHECK_PASSED'] = bool(entry['PASSED'])  # explizite Typumwandlung
        new_df['ORIG_SEQ'] = orig_seq
        new_df['TRIMMED_SEQ'] = entry['TRIMMED_SEQ']
        new_df['MASKED_SEQ'] = entry['MASKED_SEQ']
        new_df['COHORT'] = entry['COHORT']
        new_df['SUBJECT'] = entry['SUBJECT']
        new_df['TIME_POINT'] = entry['TIME_POINT']
        new_df['PLATE'] = entry['PLATE']
        new_df['PRIMER_SET'] = entry['PRIMER_SET']
        new_df['WELL'] = entry['WELL']
        new_df['SOURCE'] = entry['SOURCE']
        new_df['SUBSOURCE'] = entry['SUBSOURCE']

        # Berechne Sequenzsegmente anhand der originalen Sequenz
        # (Annahme: get_sequence und translate_nt_to_aa sind importiert und funktionsfähig)
        new_df['FWR1_NT'] = new_df.apply(lambda row: get_sequence(row['FWR1_START'], row['FWR1_END'], orig_seq, -1), axis=1)
        new_df['CDR1_NT'] = new_df.apply(lambda row: get_sequence(row['CDR1_START'], row['CDR1_END'], orig_seq, -1), axis=1)
        new_df['CDR1_AA'] = new_df['CDR1_NT'].apply(translate_nt_to_aa)
        new_df['FWR2_NT'] = new_df.apply(lambda row: get_sequence(row['FWR2_START'], row['FWR2_END'], orig_seq, -1), axis=1)
        new_df['CDR2_NT'] = new_df.apply(lambda row: get_sequence(row['CDR2_START'], row['CDR2_END'], orig_seq, -1), axis=1)
        new_df['CDR2_AA'] = new_df['CDR2_NT'].apply(translate_nt_to_aa)
        new_df['FWR3_NT'] = new_df.apply(lambda row: get_sequence(row['FWR3_START'], row['FWR3_END'], orig_seq, -1), axis=1)

        # FWR4: Falls Korrektur erforderlich, berechne FWR4_NT und FWR4_AA
        if new_df['FWR4_CORRECTION'].iloc[0]:
            new_df['FWR4_NT'] = new_df.apply(
                lambda row: get_sequence(row['FWR4_START'], row['FWR4_END'] + row['FWR4_CORRECTION'], orig_seq, -1),
                axis=1
            )
            new_df['FWR4_AA'] = new_df['FWR4_NT'].apply(translate_nt_to_aa)
        else:
            new_df['FWR4_NT'] = None
            new_df['FWR4_AA'] = None

        # Längenangaben für CDR3
        new_df['CDR3_NT'] = igblast_result.get('CDR3_NT', "N/A")
        new_df['CDR3_NT_LENGTH'] = len(cdr3_nt) if cdr3_nt != "N/A" else "N/A"
        new_df['CDR3_AA'] = igblast_result.get('CDR3_AA', "N/A")
        new_df['CDR3_AA_LENGTH'] = len(cdr3_aa) if cdr3_aa != "N/A" else "N/A"

        # Berechne C_NT und C_AA (hier wird angenommen, dass C_START im IgBLAST-Ergebnis vorliegt)
        new_df['C_NT'] = new_df.apply(lambda row: get_sequence(row['C_START'], len(orig_seq) - 20, orig_seq, -1), axis=1)
        new_df['C_AA'] = new_df['C_NT'].apply(translate_nt_to_aa)

        # Sammle diesen DataFrame
        df_list.append(new_df)

    # Führe einmaliges Concatenaten aller DataFrames durch
    combined_dataframe = pd.concat(df_list, ignore_index=True)

    # Isotyp-Bestimmung
    if isotype_determination:
        print(log_message("Determining Isotypes"), flush=True)
        combined_dataframe['ISOTYPE'] = get_isotype_with_blast(combined_dataframe['C_NT'], out_dir)
        combined_dataframe['TOP_ISOTYPE'] = combined_dataframe['ISOTYPE'].apply(lambda x: x.split("*")[0])
    else:
        combined_dataframe['ISOTYPE'] = "N.D."
        combined_dataframe['TOP_ISOTYPE'] = "N.D."

        # Falls gewünscht: Korrigiere falsche Isotyp-Spalte für Light Chains
    change_index = combined_dataframe[combined_dataframe['CHAIN_PCR'] != "HC"].index
    combined_dataframe.loc[change_index, "ISOTPYE"] = np.nan

    # Erzeuge eine B_CELL_ID anhand mehrerer Spalten
    combined_dataframe['B_CELL_ID'] = combined_dataframe[['COHORT', 'SUBJECT', 'TIME_POINT', 'TISSUE', 'SUBSET', 'PLATE', 'WELL']].apply(
        lambda x: "_".join(x.astype(str)), axis=1
    )

    # Change column order #'SELECT_HC', 'SELECT_KC', 'SELECT_LC',
    features = ['B_CELL_ID', 'SAMPLE_NAME', 'COHORT', 'SUBJECT', 'TIME_POINT',
                'TISSUE', 'SUBSET', 'PLATE', 'WELL', 'PRIMER_SET', 'SOURCE', 'SUBSOURCE', 'IGBLAST_PASSED', 'V_GENE',
                'TOP_V', 'D_GENE', 'TOP_D', 'J_GENE', 'TOP_J', 'V_IDENTITY', 'CHAIN_PCR', 'STOP_CODON', 'FRAME',
                'PRODUCTIVE', 'ORIENTATION', 'FULL_V_ALIGNMENT', 'FWR1_START', 'FWR1_END', 'FWR1_NT', 'CDR1_START',
                'CDR1_END', 'CDR1_NT', 'FWR2_START', 'FWR2_END', 'FWR2_NT', 'CDR2_START', 'CDR2_END', 'CDR2_NT', 'FWR3_START', 'FWR3_END',
                'FWR3_NT', 'CDR3_NT', 'CDR3_NT_LENGTH', 'CDR3_AA', 'CDR3_AA_LENGTH', 'V_BTOP', 'J_START', 'J_END',
                'J_BTOP', 'FWR4_START', 'FWR4_END', 'FWR4_CORRECTION', 'FWR4_NT', 'FWR4_AA', 'FWR4_FRAME',
                'FWR4_STOP_CODON', 'C_START', 'C_NT', 'C_AA', 'ISOTYPE', 'TOP_ISOTYPE', 'V_LENGTH', 'LENGTH_PASSED', 'INNER_N',
                'INNER_N_PASSED', 'MEAN_PHRED_PASSED', 'QCHECK_PASSED', 'TRIMMED_SEQ', 'MASKED_SEQ', 'ORIG_SEQ']

    combined_dataframe = combined_dataframe.reindex(columns=features)

    # Exportiere den kombinierten DataFrame als Excel-Datei
    path_to_combined_excel = Path(out_dir).joinpath(combi_ex_name)
    with pd.ExcelWriter(path_to_combined_excel) as writer:
        combined_dataframe.to_excel(writer, sheet_name="Full_Information", index=False)

    print(log_message(f"{combi_ex_name} exported to {out_dir}"), flush=True)
    return combined_dataframe


def dict_to_str(d, indent=0):
    """
    Recursively converts a dictionary into a string with indentation.

    :param d: The dictionary to convert.
    :param indent: Current indentation level (number of spaces).
    :return: A list of strings representing the dictionary.
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
    Creates a clustering log file containing a timestamp, active settings, input file names,
    subgroups to be exported, and exported file names.

    Nested dictionaries in settings are handled recursively.

    :param log_file_path: Path (string) to the log file to be created.
    :param settings: Dictionary with active settings (can be nested).
    :param input_files: List of input file names (strings).
    :param output_files: List of output file names (strings).
    :param subgroups: List of subgroups to be exported (e.g., tuples or strings).
    """
    # Create a timestamp string in the format YYMMDD-HHMMSS
    timestamp = datetime.datetime.now().strftime("%y%m%d-%H%M%S")

    with open(log_file_path, "w") as log_file:
        log_file.write(f"TIMESTAMP: {timestamp}\n\n")

        # Input files
        log_file.write("Input files:\n")
        if input_files:
            for infile in input_files:
                log_file.write(f"{infile}\n")
        else:
            log_file.write(f"- No input files passed - \n")
        log_file.write("\n")

        # Settings
        log_file.write("Active filter:\n")
        if settings:
            settings_lines = dict_to_str(settings, indent=0)
            log_file.write("\n".join(settings_lines))
        else:
            log_file.write("- No settings passed -")
        log_file.write("\n\n")

        # Output files
        log_file.write("Exported files:\n")
        if output_files:
            for outfile in output_files:
                log_file.write(f"{outfile}\n")
        else:
            log_file.write(f"- No input files passed - \n")