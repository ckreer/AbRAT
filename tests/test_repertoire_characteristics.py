import pandas as pd
import numpy as np
import pytest
from abrat.core.repertoire_characteristics import (
    get_sample_information,
    get_hydrophobicity_values,
    compute_gravy_scores,
    compute_histogram,
    compute_diversity_index,
    subsample_series_dict,
    compute_diversity_indices_for_subsampling,
    compute_pairwise_jaccard,
    sequence_similarity,
    compute_pairwise_weighted_jaccard
)

def make_sample_df():
    # Create a simple DataFrame with MultiIndex columns as expected by get_sample_information
    arrays = [
        ('SAMPLE_INFORMATION', 'COHORT'),
        ('SAMPLE_INFORMATION', 'SUBJECT'),
        ('SAMPLE_INFORMATION', 'TIME_POINT'),
        ('SAMPLE_INFORMATION', 'TISSUE'),
        ('SAMPLE_INFORMATION', 'SUBSET'),
        ('HEAVY_CHAIN', 'HEAVY_FOUND'),
        ('LIGHT_CHAIN', 'LIGHT_FOUND'),
        ('LIGHT_CHAIN', 'PCR_ISOTYPE')
    ]
    tuples = arrays
    idx = pd.MultiIndex.from_tuples(tuples)
    data = [
        ['C1', 'D1', 'T1', 'S1', 'SS1', True, True, 'KC'],
        ['C1', 'D1', 'T1', 'S1', 'SS1', False, True, 'LC']
    ]
    df = pd.DataFrame(data, columns=idx)
    return df

def test_get_sample_information():
    df = make_sample_df()
    info = get_sample_information(df)
    assert info['Cohort'] == 'C1'
    assert info['Donor'] == 'D1'
    assert info['Time points'] == 'T1'
    assert info['Sample'] == 'S1'
    assert info['Subset'] == 'SS1'
    assert info['Total BCRs'] == 2
    assert info['Heavy chains'] == 1
    assert info['Light chains'] == 2
    assert info['Kappa chains'] == 1
    assert info['Lambda chains'] == 1

def test_get_hydrophobicity_values_and_error():
    seq = "ACDE"
    vals = get_hydrophobicity_values(seq, scale='kyte-doolittle')
    # Known values length
    assert len(vals) == 4
    with pytest.raises(ValueError):
        get_hydrophobicity_values("ACBX", scale='kyte-doolittle')
    with pytest.raises(ValueError):
        get_hydrophobicity_values("ACDE", scale='unknown')

def test_compute_gravy_scores():
    series = pd.Series(["ACD", "GGG"])
    gravy = compute_gravy_scores(series, scale='eisenberg')
    # 'ACD' -> [0.62, 0.29, -0.90]
    assert isinstance(gravy, pd.Series)
    assert np.isclose(gravy.iloc[0], np.mean([0.62, 0.29, -0.90]))
    # 'GGG' -> three glycine values (0.48 each)
    assert np.isclose(gravy.iloc[1], 0.48)

def test_compute_histogram_and_errors():
    series = pd.Series([0.1, 0.5, 1.2, 1.8, 2.0])
    hist = compute_histogram(series, bin_size=1)
    # bins [0.1,1.1), [1.1,2.1)
    assert hist.sum() == 5
    with pytest.raises(ValueError):
        compute_histogram(pd.Series(['a','b']), bin_size=1)

def test_compute_diversity_index():
    seqs = pd.Series(['A','A','B','C'])
    # Shannon: -[0.5*ln0.5 + 0.25*ln0.25 + 0.25*ln0.25]
    shannon = compute_diversity_index(seqs, index_type='shannon')
    assert pytest.approx(shannon) == - (0.5*np.log(0.5) + 0.25*np.log(0.25) + 0.25*np.log(0.25))
    inv_simpson = compute_diversity_index(seqs, index_type='inverse_simpson')
    assert np.isclose(inv_simpson, 1 / (0.5**2 + 0.25**2 + 0.25**2))
    with pytest.raises(ValueError):
        compute_diversity_index(seqs, index_type='unknown')

def test_subsample_and_diversity_indices_for_subsampling():
    d = {'a': pd.Series(list('AAA')), 'b': pd.Series(list('ABC'))}
    subs, min_size, key = subsample_series_dict(d)
    assert min_size == 3
    # generate diversity dict
    res = compute_diversity_indices_for_subsampling(d, index_type='shannon', iterations=5, random_state=0)
    assert set(res.keys()) == {'a','b'}
    assert all('mean' in res[k] and 'std' in res[k] for k in res)

def test_pairwise_jaccard_and_weighted():
    d = {'x': ['A','B','C'], 'y': ['B','C','D']}
    results, matrix = compute_pairwise_jaccard(d, iterations=5, random_state=0)
    assert ('x','y') in results
    # Weighted Jaccard
    wres, wmat = compute_pairwise_weighted_jaccard(d, iterations=5, random_state=0)
    assert ('x','y') in wres
    assert np.allclose(matrix.loc['x','y'], results[('x','y')]['mean'])
    assert isinstance(wmat, pd.DataFrame)

def test_sequence_similarity():
    assert sequence_similarity("ABC","ABC") == 1.0
    assert 0 <= sequence_similarity("ABC","DEF") <= 1.0
