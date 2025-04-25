import os
import pytest
from pathlib import Path

import abrat.core.build_igblast_db as module

def test_is_segment_functions():
    assert module.is_v_segment(">IGHV1-18 description")
    assert module.is_v_segment(">IGKV3-20")
    assert module.is_v_segment(">IGLV1-2")
    assert not module.is_v_segment(">IGHD123")
    assert module.is_d_segment(">IGHD123")
    assert not module.is_d_segment(">IGHV123")
    assert module.is_j_segment(">IGHJ1")
    assert module.is_j_segment(">IGKJ2")
    assert module.is_j_segment(">IGLJ3")
    assert not module.is_j_segment(">IGHD4")

def test_split_fasta_by_segment(tmp_path):
    # create a FASTA file with one of each type
    fasta = tmp_path / "test.fasta"
    fasta.write_text(
        ">IGHV1\nAAA\n"
        ">IGHD1\nTTT\n"
        ">IGHJ1\nGGG\n"
    )
    v = tmp_path / "v.fasta"
    d = tmp_path / "d.fasta"
    j = tmp_path / "j.fasta"
    module.split_fasta_by_segment([str(fasta)], str(v), str(d), str(j))
    assert ">IGHV1" in v.read_text()
    assert ">IGHD1" in d.read_text()
    assert ">IGHJ1" in j.read_text()

def test_download_file(tmp_path, monkeypatch):
    # simulate requests.get
    class DummyResponse:
        def raise_for_status(self): pass
        def iter_content(self, chunk_size):
            yield b'foo'
            yield b'bar'
    # Monkey-patch module.requests.get
    monkeypatch.setattr(module, 'requests', type('R', (), {'get': lambda *args, **kwargs: DummyResponse()}))
    target = tmp_path / "out.bin"
    module.download_file("http://example.com/file", str(target))
    assert target.read_bytes() == b'foobar'

def test_make_blast_db(monkeypatch):
    calls = []
    def dummy_run(cmd, check):
        calls.append(cmd)
    # Monkey-patch subprocess.run
    monkeypatch.setattr(module.subprocess, 'run', dummy_run)
    module.make_blast_db("in.fasta", "dbout")
    assert calls, "subprocess.run was not called"
    cmd = calls[0]
    assert cmd[0] == 'makeblastdb'
    assert '-in' in cmd and 'in.fasta' in cmd
    assert '-out' in cmd and 'dbout' in cmd
