import pandas as pd
import numpy as np
import pytest
from abrat.core.clonal_assignment import (
    parse_clone_number,
    normalized_levenshtein,
    standardize_missing_values,
    iterative_greedy_cdr3_clustering,
    matrix_greedy_cdr3_clustering,
    hierarchical_cdr3_clustering,
)

def test_parse_clone_number():
    assert parse_clone_number("Clone-1") == 1
    assert parse_clone_number("Clone-123") == 123
    assert parse_clone_number("Foo-1") == float('inf')

def test_normalized_levenshtein():
    assert normalized_levenshtein("", "") == 0
    assert normalized_levenshtein("A", "A") == pytest.approx(0.0)
    assert normalized_levenshtein("A", "T") == pytest.approx(1.0)
    assert normalized_levenshtein("AA", "AT") == pytest.approx(0.5)

def test_standardize_missing_values():
    df = pd.DataFrame({'col': ["N/A", "nan", "", "valid"]})
    out = standardize_missing_values(df)
    assert pd.isna(out.loc[0, 'col'])
    assert pd.isna(out.loc[1, 'col'])
    assert pd.isna(out.loc[2, 'col'])
    assert out.loc[3, 'col'] == "valid"

def test_iterative_greedy_basic():
    df = pd.DataFrame({'CDR3_AA': ['AAA', 'AAA', 'TTT', None]}, index=['a','b','c','d'])
    clusters, reps = iterative_greedy_cdr3_clustering(df, {'length_threshold':0, 'lev_threshold':0.0})
    # expect two non-zero clusters: cluster 1 for 'AAA', cluster 2 for 'TTT'
    # cluster labels start at 1 for the first cluster
    assert set(clusters[clusters>0]) == {1, 2}
    assert reps['a'] in ['a', 'b']
    assert clusters['d'] == 0

def test_matrix_greedy_basic():
    df = pd.DataFrame({'CDR3_AA': ['AAA', 'AAA', 'TTT', None]}, index=['a','b','c','d'])
    clusters, reps = matrix_greedy_cdr3_clustering(df, {'length_threshold':0, 'lev_threshold':0.0})
    # expected clusters: 0 for 'AAA', 1 for 'TTT'
    assert set(clusters.unique()) == {0, 1}
    assert set(clusters[clusters>0]) == {1}
    assert reps['c'] == 'c'
    assert clusters['d'] == 0

def test_hierarchical_single_and_empty():
    # single valid
    df1 = pd.DataFrame({'CDR3_AA': ['AAA']}, index=['x'])
    cl1, rp1 = hierarchical_cdr3_clustering(df1, {'lev_threshold':0.5, 'length_threshold':1})
    assert cl1['x'] == 1
    assert rp1['x'] == 'x'
    # empty
    df2 = pd.DataFrame({'CDR3_AA': []})
    cl2, rp2 = hierarchical_cdr3_clustering(df2, {})
    assert cl2.empty
    assert rp2.empty
