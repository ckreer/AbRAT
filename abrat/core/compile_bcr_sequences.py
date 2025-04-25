#!/usr/bin/env python3

import pandas as pd
import numpy as np

def edit_distance(s1, s2, exclude_list):
    """
    Computes a Hamming-like edit distance between two strings, ignoring positions where either letter is in the exclude list.

    Parameters:
        s1 (str): The first string.
        s2 (str): The second string.
        exclude_list (list): List of characters to exclude from the comparison.

    Returns:
        int: The computed edit distance.
    """
    return sum(
        [pair[0] != pair[1] for pair in list(zip(s1, s2))
         if not any([element in pair for element in exclude_list])]
    )


def find_best_sequence(df, difference_cut_off):
    """
    Determines the best sequence hit from a DataFrame based on the 'INNER_N' metric and performs pairwise comparisons.

    The function selects the best hit by sorting the DataFrame by the 'INNER_N' column.
    If multiple entries share the minimal 'INNER_N', it prefers the "MINI" source, and if necessary,
    then the "MIDI" source. It also performs pairwise comparisons (using masked sequences, V gene, and J gene)
    among quality-checked entries to identify collisions when differences exceed the specified cutoff.

    Parameters:
        df (pd.DataFrame): DataFrame containing sequence data.
        difference_cut_off (int): Cutoff value for differences to identify collisions.

    Returns:
        tuple: A tuple (best_df, sources, others, collisions) where:
            - best_df (pd.DataFrame): DataFrame with the best sequence hit.
            - sources (list): List of unique sources found in the DataFrame.
            - others (str): Comma-separated string of sample names excluding the best hit.
            - collisions (list): List of tuples describing collisions between sequences.
    """
    # TODO: How to proceed with truncated sequences?
    df = df[df['INNER_N'] != 'N/A']
    sources = list(df['SOURCE'].unique())
    collisions = []
    best_df = df.sort_values(by=['INNER_N']).head(1)
    if len(df) > 1:
        # Find the minimal INNER_N and isolate all entries with that minimal value.
        min_n_df = df[df['INNER_N'] == df['INNER_N'].min()]
        # If several hits with the same INNER_N exist, prefer "MINI" over others.
        if (len(min_n_df) > 1) and ("MINI" in list(min_n_df['SOURCE'].unique())):
            best_df = min_n_df[min_n_df['SOURCE'] == "MINI"].head(1)
        if (len(min_n_df) > 1) and ("MIDI" in list(min_n_df['SOURCE'].unique())):
            best_df = min_n_df[min_n_df['SOURCE'] == "MIDI"].head(1)
        best = best_df['SAMPLE_NAME'].values[0]
        name_list = list(df['SAMPLE_NAME'])
        name_list.remove(best)
        others = ", ".join(name_list)
        # Pairwise comparison of sequences, V gene, and J gene.
        for idx1 in df[df['QCHECK_PASSED']].index[:-1]:
            seq1 = df.at[idx1, 'MASKED_SEQ']
            vg1 = df.at[idx1, 'V_GENE'].split("*")[0]
            jg1 = df.at[idx1, 'J_GENE'].split("*")[0]
            for idx2 in df[df['QCHECK_PASSED']].index[1:]:
                seq2 = df.at[idx2, 'MASKED_SEQ']
                vg2 = df.at[idx2, 'V_GENE'].split("*")[0]
                jg2 = df.at[idx2, 'J_GENE'].split("*")[0]
                vgs = "IDENTICAL_V_GENE" if vg1 == vg2 else "DIFFERENT_V_GENE"
                jgs = "IDENTICAL_J_GENE" if jg1 == jg2 else "DIFFERENT_J_GENE"
                distance = edit_distance(seq1, seq2, ["N"])
                if (distance > difference_cut_off) or (vgs == "DIFFERENT_V_GENE") or (jgs == "DIFFERENT_J_GENE"):
                    collisions.append((df.at[idx1, 'SAMPLE_NAME'], df.at[idx2, 'SAMPLE_NAME'], vgs, jgs,
                                       "DISTANCE(NT)=" + str(distance)))
    else:
        others = np.nan
    return best_df, sources, others, collisions


def add_to_list(l, value):
    """
    Adds a value to a list and returns the updated list. If the provided list is None or evaluates to False,
    a new list is created.

    Parameters:
        l (list or pd.Series): The original list or series.
        value: The value to add.

    Returns:
        list: The updated list with the value appended.
    """
    if l.iloc[0]:
        l += [", " + str(value)]
    else:
        l = [value]
    return l


def compile_bcrs(composite_df, h_ident, k_ident, l_ident, col_co):
    """
    Extracts heavy (H), kappa (K), and lambda (L) chain information per identifier and returns a combined DataFrame
    with the best sequence selected according to the find_best_sequence function.

    The function groups the composite DataFrame by 'B_CELL_ID'. For each group, it retrieves all heavy, kappa, and lambda
    chain sequences, determines the best sequence for each chain using find_best_sequence, assigns warnings for collisions
    and other issues, and consolidates all relevant chain information into a single DataFrame.

    Parameters:
        composite_df (pd.DataFrame): Input DataFrame containing combined BCR information.
        h_ident (str): Identifier for heavy chains.
        k_ident (str): Identifier for kappa chains.
        l_ident (str): Identifier for lambda chains.
        col_co (int): Cutoff value for the 'INNER_N' metric used in find_best_sequence.

    Returns:
        pd.DataFrame: A DataFrame with the compiled antibody (BCR) information.
    """
    if composite_df.empty:
        raise ValueError("Cannot compile BCRs from empty dataframe.")
    sample_columns = ['B_CELL_ID', 'COHORT', 'SUBJECT', 'TIME_POINT', 'TISSUE',
                      'SUBSET', 'PLATE', 'WELL']
    clone_columns = ['HC-LC_CLUSTER', 'CLUSTER_SIZE', 'IS_CLONAL', 'CLONE', 'CLONE_COLOR']

    hc_columns = ['HEAVY_FOUND', 'ISOTYPE', 'TOP_ISOTYPE', 'HC_SUBCLUSTER', 'CLUSTER_REPRESENTATIVE',
                  'SAMPLE_NAME', 'SOURCE', 'SUBSOURCE', 'V_GENE', 'TOP_V', 'D_GENE', 'TOP_D', 'J_GENE', 'TOP_J',
                  'STOP_CODON', 'FRAME', 'PRODUCTIVE', 'FULL_V_ALIGNMENT', 'CDR3_AA', 'CDR3_AA_LENGTH',
                  'CDR3_NT', 'CDR3_NT_LENGTH', 'V_IDENTITY', 'V_BTOP', 'INNER_N', 'QCHECK_PASSED', 'MIDI_WARNING',
                  'COLLISIONS', 'FWR4_WARNING', 'ALTERNATIVE_SEQ', 'ORIGINAL_SEQ', 'TRIMMED_SEQ', 'MASKED_SEQ']
    light_chain_columns = ['LIGHT_FOUND', 'SAMPLE_NAME', 'SOURCE', 'SUBSOURCE', 'V_GENE', 'TOP_V',
                           'J_GENE', 'TOP_J', 'STOP_CODON', 'FRAME', 'PRODUCTIVE', 'FULL_V_ALIGNMENT',
                           'CDR3_AA', 'CDR3_AA_LENGTH', 'CDR3_NT', 'CDR3_NT_LENGTH', 'V_IDENTITY', 'V_BTOP',
                           'INNER_N', 'QCHECK_PASSED', 'MIDI_WARNING', 'COLLISIONS', 'FWR4_WARNING',
                           'ALTERNATIVE_SEQ', 'ORIGINAL_SEQ', 'TRIMMED_SEQ', 'MASKED_SEQ']
    kc_columns = ['KAPPA_FOUND', 'SAMPLE_NAME', 'SOURCE', 'SUBSOURCE', 'V_GENE', 'J_GENE',
                  'STOP_CODON', 'FRAME', 'PRODUCTIVE', 'FULL_V_ALIGNMENT', 'CDR3_AA', 'CDR3_AA_LENGTH',
                  'CDR3_NT', 'CDR3_NT_LENGTH', 'V_IDENTITY', 'V_BTOP', 'INNER_N', 'QCHECK_PASSED', 'MIDI_WARNING',
                  'COLLISIONS', 'FWR4_WARNING', 'ALTERNATIVE_SEQ', 'ORIGINAL_SEQ', 'TRIMMED_SEQ', 'MASKED_SEQ']
    lc_columns = ['LAMBDA_FOUND', 'SAMPLE_NAME', 'SOURCE', 'SUBSOURCE', 'V_GENE', 'J_GENE', 'V_BTOP',
                  'STOP_CODON', 'FRAME', 'PRODUCTIVE', 'FULL_V_ALIGNMENT', 'CDR3_AA', 'CDR3_AA_LENGTH',
                  'CDR3_NT', 'CDR3_NT_LENGTH', 'V_IDENTITY', 'INNER_N', 'QCHECK_PASSED', 'MIDI_WARNING',
                  'COLLISIONS', 'FWR4_WARNING', 'ALTERNATIVE_SEQ', 'ORIGINAL_SEQ', 'TRIMMED_SEQ', 'MASKED_SEQ']

    fin_df = pd.DataFrame()

    # For each B_CELL_ID, retrieve all heavy, kappa, and lambda chain sequences, determine the best sequence, and save alternative sequences.
    for individual_b_cell, individual_b_cell_df in composite_df.groupby('B_CELL_ID'):
        # For each well, find all heavy, kappa, and lambda sequences, determine the best sequence, and record alternative sequences.
        sample_infos = individual_b_cell_df[sample_columns].head(1).set_index('B_CELL_ID')
        sample_infos["WARNING"] = ""

        best_h, h_sources, other_hs, h_collisions = find_best_sequence(
            individual_b_cell_df[individual_b_cell_df['CHAIN_PCR'] == h_ident], col_co)
        best_k, k_sources, other_ks, k_collisions = find_best_sequence(
            individual_b_cell_df[individual_b_cell_df['CHAIN_PCR'] == k_ident], col_co)
        best_l, l_sources, other_ls, l_collisions = find_best_sequence(
            individual_b_cell_df[individual_b_cell_df['CHAIN_PCR'] == l_ident], col_co)

        h_found = False if best_h.empty else True
        k_found = False if best_k.empty else True
        l_found = False if best_l.empty else True

        best_h['ALTERNATIVE_SEQ'] = other_hs
        best_k['ALTERNATIVE_SEQ'] = other_ks
        best_l['ALTERNATIVE_SEQ'] = other_ls

        # Add collision warnings to sample info.
        if h_collisions:
            best_h['COLLISIONS'] = [h_collisions]
            sample_infos['WARNING'] = add_to_list(sample_infos['WARNING'], 'COLLISION_HC')
        else:
            best_h['COLLISIONS'] = None

        if k_collisions:
            best_k['COLLISIONS'] = [k_collisions]
            sample_infos['WARNING'] = add_to_list(sample_infos['WARNING'], 'COLLISION_KC')
        else:
            best_k['COLLISIONS'] = None

        if l_collisions:
            best_l['COLLISIONS'] = [l_collisions]
            sample_infos['WARNING'] = add_to_list(sample_infos['WARNING'], 'COLLISION_LC')
        else:
            best_l['COLLISIONS'] = None

        # Add MIDI warnings if the MIDI source is present but not selected as the best sequence.
        if ("MIDI" in h_sources) and not ("MIDI" in best_h['SOURCE'].values):
            best_h['MIDI_WARNING'] = "MIDI not the best Sequence"
            sample_infos['WARNING'] = add_to_list(sample_infos['WARNING'], 'MIDI_WARNING_HC')
        else:
            best_h['MIDI_WARNING'] = None

        if ("MIDI" in k_sources) and not ("MIDI" in best_k['SOURCE'].values):
            best_k['MIDI_WARNING'] = "MIDI not the best Sequence"
            sample_infos['WARNING'] = add_to_list(sample_infos['WARNING'], 'MIDI_WARNING_KC')
        else:
            best_k['MIDI_WARNING'] = None

        if ("MIDI" in l_sources) and not ("MIDI" in best_l['SOURCE'].values):
            best_l['MIDI_WARNING'] = "MIDI not the best Sequence"
            sample_infos['WARNING'] = add_to_list(sample_infos['WARNING'], 'MIDI_WARNING_LC')
        else:
            best_l['MIDI_WARNING'] = None

        # Add warnings if FWR4 is out-of-frame or contains a stop codon.
        for chain, best_chain_df in (('HC', best_h), ('KC', best_k), ('LC', best_l)):
            chain_fwr4_warnings = []
            if "Out-of-frame" in best_chain_df['FWR4_FRAME']:
                chain_fwr4_warnings.append('FWR4 is out-of-frame')
            if "Yes" in best_chain_df['FWR4_STOP_CODON']:
                chain_fwr4_warnings.append('FWR4 stop codon found')
            if not chain_fwr4_warnings:
                chain_fwr4_warnings = None
            if chain_fwr4_warnings:
                sample_infos['WARNING'] = add_to_list(sample_infos['WARNING'], 'FWR4_WARNING_' + chain)
            best_chain_df['FWR4_WARNING'] = None if chain_fwr4_warnings is None else ", ".join(chain_fwr4_warnings)

        # Determine the best light chain: if both kappa and lambda are found, select based on productivity and minimal INNER_N.
        best_light_chain = pd.DataFrame(columns=light_chain_columns)
        best_light_chain['LIGHT_FOUND'] = False
        if k_found and not l_found:
            best_light_chain = best_k.reindex(light_chain_columns, axis=1)
            best_light_chain['LIGHT_FOUND'] = True
            best_light_chain['PCR_ISOTYPE'] = "KC"
        elif l_found and not k_found:
            best_light_chain = best_l.reindex(light_chain_columns, axis=1)
            best_light_chain['LIGHT_FOUND'] = True
            best_light_chain['PCR_ISOTYPE'] = "LC"
        elif not k_found and not l_found:
            sample_infos['WARNING'] = add_to_list(sample_infos['WARNING'], 'NO_LIGHT_CHAINS')
            best_light_chain['LIGHT_FOUND'] = False
        elif k_found and l_found:
            best_light_chain['LIGHT_FOUND'] = True
            prod_k = True if str(best_k['PRODUCTIVE'].values[0]) == "Yes" else False
            prod_l = True if str(best_l['PRODUCTIVE'].values[0]) == "Yes" else False
            # If both are productive.
            if prod_k and prod_l:
                sample_infos['WARNING'] = add_to_list(sample_infos['WARNING'], 'TWO_PRODUCTIVE_LIGHT_CHAINS')
                qcheck_k = False if type(best_k['QCHECK_PASSED'].values[0]) == str else best_k['QCHECK_PASSED'].values[0]
                qcheck_l = False if type(best_l['QCHECK_PASSED'].values[0]) == str else best_l['QCHECK_PASSED'].values[0]
                if qcheck_k and not qcheck_l:
                    best_light_chain = best_k.reindex(light_chain_columns, axis=1)
                    best_light_chain['PCR_ISOTYPE'] = "KC"
                elif not qcheck_k and qcheck_l:
                    best_light_chain = best_l.reindex(light_chain_columns, axis=1)
                    best_light_chain['PCR_ISOTYPE'] = "LC"
                else:
                    # Choose the sequence with fewer 'N's if both or none passed quality check.
                    if int(best_k['INNER_N'].fillna(1000).iloc[0]) <= int(best_l['INNER_N'].fillna(1000).iloc[0]):
                        best_light_chain = best_k.reindex(light_chain_columns, axis=1)
                        best_light_chain['PCR_ISOTYPE'] = "KC"
                    else:
                        best_light_chain = best_l.reindex(light_chain_columns, axis=1)
                        best_light_chain['PCR_ISOTYPE'] = "LC"
            # If only one is productive.
            elif prod_k and not prod_l:
                best_light_chain = best_k.reindex(light_chain_columns, axis=1)
                best_light_chain['PCR_ISOTYPE'] = "KC"
                sample_infos['WARNING'] = add_to_list(sample_infos['WARNING'], 'TWO_LIGHT_CHAINS_ONE_PRODUCTIVE')
            elif not prod_k and prod_l:
                best_light_chain = best_l.reindex(light_chain_columns, axis=1)
                best_light_chain['PCR_ISOTYPE'] = "LC"
                sample_infos['WARNING'] = add_to_list(sample_infos['WARNING'], 'TWO_LIGHT_CHAINS_ONE_PRODUCTIVE')
            # If both are unproductive, select the one with fewer 'N's.
            elif not prod_k and not prod_l:
                if int(best_k['INNER_N'].fillna(1000).iloc[0]) <= int(best_l['INNER_N'].fillna(1000).iloc[0]):
                    best_light_chain = best_k.reindex(light_chain_columns, axis=1)
                    best_light_chain['PCR_ISOTYPE'] = "KC"
                else:
                    best_light_chain = best_l.reindex(light_chain_columns, axis=1)
                    best_light_chain['PCR_ISOTYPE'] = "LC"
                sample_infos['WARNING'] = add_to_list(sample_infos['WARNING'], 'TWO_UNPRODUCTIVE_LIGHT_CHAINS')

        # Re-index the best chain DataFrames.
        best_h = best_h.set_index('B_CELL_ID')
        best_k = best_k.set_index('B_CELL_ID')
        best_l = best_l.set_index('B_CELL_ID')
        best_light_chain['B_CELL_ID'] = individual_b_cell
        best_light_chain = best_light_chain.set_index('B_CELL_ID')

        best_h['HEAVY_FOUND'] = False
        best_k['KAPPA_FOUND'] = False
        best_l['LAMBDA_FOUND'] = False

        best_h.at[individual_b_cell, 'HEAVY_FOUND'] = h_found  # Changed on 2024-11-28; originally encoded as "Yes"/"No"
        best_light_chain['LIGHT_FOUND'] = any([k_found, l_found])  # Changed on 2024-11-28; originally encoded as "Yes"/"No"
        best_k.at[individual_b_cell, 'KAPPA_FOUND'] = k_found  # Changed on 2024-11-28; originally encoded as "Yes"/"No"
        best_l.at[individual_b_cell, 'LAMBDA_FOUND'] = l_found  # Changed on 2024-11-28; originally encoded as "Yes"/"No"

        sample_infos["CHAINS"] = ", ".join(
            [["H", "K", "L"][i] for i, chain in enumerate([h_found, k_found, l_found]) if chain]
        ) if any([h_found, k_found, l_found]) else "NO CHAINS FOUND"

        # Add blank columns for clonal analysis.
        for col in clone_columns:
            sample_infos[col] = pd.Series()

        best_h['HC_SUBCLUSTER'] = pd.Series()
        best_h['CLUSTER_REPRESENTATIVE'] = pd.Series()
        best_light_chain['LC_SUBCLUSTER'] = pd.Series()
        best_light_chain['CLUSTER_REPRESENTATIVE'] = pd.Series()

        # Add blank columns for chain selection (for subsequent cloning).
        new_columns = ['SELECT_HC', 'SELECT_KC', 'SELECT_LC']
        sample_infos.loc[:, new_columns] = np.nan

        # Change column order: place the index column first, then new columns, followed by remaining columns.
        first_col = [sample_infos.columns[0]]
        remaining_columns = [col for col in sample_infos.columns if col not in set(first_col).union(new_columns)]
        new_order = first_col + new_columns + remaining_columns

        # Concatenate all required information.
        antibody_df = pd.concat(
            [sample_infos[new_order], best_h[hc_columns],
             best_light_chain[light_chain_columns],
             best_k[kc_columns], best_l[lc_columns]],
            axis=1,
            keys=['SAMPLE_INFORMATION', 'HEAVY_CHAIN', 'LIGHT_CHAIN', 'KAPPA_CHAIN', 'LAMBDA_CHAIN']
        )
        # Append to the final DataFrame.
        fin_df = pd.concat([fin_df, antibody_df], axis=0)

    return fin_df
