import pandas as pd
import numpy as np
import pytest
from abrat.core.compile_bcr_sequences import (
    edit_distance,
    find_best_sequence,
    add_to_list,
    compile_bcrs
)

def test_edit_distance():
    assert edit_distance("ABC", "ABC", []) == 0
    assert edit_distance("ABC", "ABD", []) == 1
    # Excluding 'N' so mismatch at position with 'N' ignored
    assert edit_distance("ANC", "ABC", ["N"]) == 0

def test_find_best_sequence_single():
    df = pd.DataFrame([{
        'INNER_N': 1,
        'SOURCE': 'X',
        'SAMPLE_NAME': 's1',
        'QCHECK_PASSED': True,
        'MASKED_SEQ': 'AAA',
        'V_GENE': 'V*01',
        'J_GENE': 'J*01'
    }])
    best_df, sources, others, collisions = find_best_sequence(df, difference_cut_off=0)
    # Single entry: best_df is same df
    assert best_df.iloc[0]['SAMPLE_NAME'] == 's1'
    assert sources == ['X']
    assert pd.isna(others)
    assert collisions == []

def test_find_best_sequence_multiple_no_collision():
    df = pd.DataFrame([
        {'INNER_N': 1, 'SOURCE': 'A', 'SAMPLE_NAME': 's1', 'QCHECK_PASSED': True,
         'MASKED_SEQ': 'AAA', 'V_GENE': 'V*01', 'J_GENE': 'J*01'},
        {'INNER_N': 2, 'SOURCE': 'B', 'SAMPLE_NAME': 's2', 'QCHECK_PASSED': True,
         'MASKED_SEQ': 'AAA', 'V_GENE': 'V*01', 'J_GENE': 'J*01'},
    ])
    best_df, sources, others, collisions = find_best_sequence(df, difference_cut_off=0)
    # Best is s1
    assert best_df.iloc[0]['SAMPLE_NAME'] == 's1'
    assert set(sources) == {'A', 'B'}
    assert others == 's2'
    # identical sequences => no collisions
    assert collisions == []

def test_add_to_list():
    # Starting from empty string in series
    l = pd.Series([""])
    out = add_to_list(l, 'X')
    assert out == ['X']

def test_compile_bcrs_empty():
    empty = pd.DataFrame()
    with pytest.raises(ValueError) as ei:
        compile_bcrs(empty, 'H', 'K', 'L', col_co=0)
    assert "Cannot compile BCRs from empty dataframe" in str(ei.value)
