import pytest
import os
from pathlib import Path
from abrat.core.analyze_ab1_files import combine_ab1_files, run_igblast, get_infos_from_igblast

def test_combine_ab1_files_empty(tmp_path):
    empty_dir = tmp_path / "empty_dir"
    empty_dir.mkdir()
    result = combine_ab1_files(empty_dir, tmp_path, "out.fasta")
    assert isinstance(result, tuple)
    ok, msg = result
    assert not ok
    assert "No .ab1-files found" in msg

def test_combine_ab1_files_nonexistent(tmp_path):
    nonexist = tmp_path / "does_not_exist"
    result = combine_ab1_files(nonexist, tmp_path, "out.fasta")
    assert isinstance(result, tuple)
    ok, msg = result
    assert not ok
    assert "Error:" in msg

def test_run_igblast_input_not_exist(tmp_path):
    out_dir = tmp_path / "out"
    ok, msg = run_igblast("nonexistent.file", out_dir, "fmt")
    assert not ok
    assert "does not exist" in msg

def test_get_infos_minimal():
    blast = "# Query: sample123\n"
    name, features = get_infos_from_igblast(blast)
    assert name == "sample123"
    assert isinstance(features, dict)
    for key in ("V_GENE", "J_GENE", "CDR3_NT", "IGBLAST_PASSED"):
        assert key in features