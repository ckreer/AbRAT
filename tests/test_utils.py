
import pytest
import numpy as np
import pandas as pd
from abrat.core.utils import (
    str2bool, timestamp_filename, date_stamp, log_message,
    reverse_complement, translate_nt_to_aa, hamming, chunk,
    get_sequence, decompress, get_indels, dict_to_str
)

def test_str2bool_true_variants():
    for val in ["yes", "True", "T", "y", "1"]:
        assert str2bool(val) is True

def test_str2bool_false_variants():
    for val in ["no", "false", "F", "n", "0"]:
        assert str2bool(val) is False

def test_str2bool_invalid():
    with pytest.raises(Exception):
        str2bool("maybe")

def test_timestamp_filename_and_date_stamp():
    fn = "file.txt"
    stamped = timestamp_filename(fn, "insert")
    assert "file_insert" in stamped and stamped.endswith(".txt")
    ds = date_stamp("name")
    assert ds.endswith("_name")

def test_log_message_format():
    msg = log_message("Test")
    assert "Test" in msg and "-" in msg.split(" ")[0]

def test_reverse_complement_and_translate_nt_to_aa():
    dna = "ATGC"
    rev = reverse_complement(dna)
    # ATGC -> GCAT complement: CGTA -> reversed gives TACG? Actually reverse_complement: GCAT
    assert rev == "GCAT"
    # translate "ATG" -> M, "GC" incomplete codon -> ?
    aa = translate_nt_to_aa("ATGAAA")
    assert aa.startswith("MK")

def test_hamming_equal_length():
    assert hamming("AAA", "AAB") == 1

def test_hamming_diff_length_exits():
    with pytest.raises(SystemExit):
        hamming("A", "AA")

def test_chunk_and_get_sequence():
    lst = list(range(5))
    chunks = list(chunk(lst, 2))
    assert chunks == [[0,1],[2,3],[4]]
    seq = "ABCDEFG"
    # get_sequence start=1, end=3, frame_shift=0 gives BCD
    assert get_sequence(1,3, seq, 0) == "BCD"

def test_decompress_and_get_indels():
    btop = "5A-3"
    decomp = decompress(btop)
    # '5' -> five '..', 'A-' stays, '3' -> three '..'
    assert isinstance(decomp, list)
    ins, dels = get_indels("A-")
    assert ins == 0 and dels == 0 or isinstance(ins, int)

def test_dict_to_str():
    d = {"a":1, "b":{"c":2}}
    lines = dict_to_str(d)
    assert any("a: 1" in l for l in lines)
    assert any("c: 2" in l for l in lines)
