"""Tests for dna_tools.py - run with: python3 test_dna_tools.py"""
import subprocess
import sys
import tempfile
from pathlib import Path

HERE = Path(__file__).parent
sys.path.insert(0, str(HERE))

import dna_tools as dt  # noqa: E402


def test_gc_content_basic_and_n_excluded():
    assert dt.gc_content("GGCC") == 100.0
    assert dt.gc_content("ATAT") == 0.0
    assert dt.gc_content("ATGC") == 50.0
    assert dt.gc_content("GCNN") == 100.0  # Ns excluded from the denominator
    print("PASS gc_content including N exclusion")


def test_reverse_complement():
    assert dt.reverse_complement("ATGC") == "GCAT"
    assert dt.reverse_complement("AAATTTNN") == "NNAAATTT"
    print("PASS reverse_complement")


def test_transcribe():
    assert dt.transcribe("ATGC") == "AUGC"
    assert dt.transcribe("TTN") == "UUN"
    print("PASS transcribe")


def test_find_motif_positions_and_overlaps():
    assert dt.find_motif("ATGATG", "ATG") == [1, 4]
    assert dt.find_motif("AAAA", "AA") == [1, 2, 3]  # overlapping hits
    assert dt.find_motif("CCCC", "ATG") == []
    print("PASS find_motif (1-based, overlapping)")


def test_read_fasta_valid_and_errors():
    with tempfile.TemporaryDirectory() as d:
        good = Path(d) / "g.fasta"
        good.write_text(">one\natgc\nNN\n>two\nTTTT\n")
        recs = dt.read_fasta(good)
        assert recs == {"one": "ATGCNN", "two": "TTTT"}  # uppercased, multi-line joined

        bad_base = Path(d) / "b.fasta"
        bad_base.write_text(">x\nACGTX\n")
        try:
            dt.read_fasta(bad_base)
            raise AssertionError("accepted invalid base")
        except ValueError as e:
            assert "invalid base" in str(e).lower()

        no_header = Path(d) / "n.fasta"
        no_header.write_text("ACGT\n")
        try:
            dt.read_fasta(no_header)
            raise AssertionError("accepted headerless file")
        except ValueError as e:
            assert "FASTA" in str(e)

        empty = Path(d) / "e.fasta"
        empty.write_text("")
        try:
            dt.read_fasta(empty)
            raise AssertionError("accepted empty file")
        except ValueError as e:
            assert "No sequences" in str(e)
    print("PASS read_fasta (valid, invalid base, missing header, empty file)")


def test_cli_end_to_end():
    demo = HERE / "demo.fasta"
    r = subprocess.run([sys.executable, str(HERE / "dna_tools.py"), "stats", str(demo)],
                       capture_output=True, text=True)
    assert r.returncode == 0 and "demo_gene_fragment" in r.stdout and "GC" in r.stdout
    r = subprocess.run([sys.executable, str(HERE / "dna_tools.py"), "motif", str(demo), "ATG"],
                       capture_output=True, text=True)
    assert r.returncode == 0 and "position(s)" in r.stdout
    r = subprocess.run([sys.executable, str(HERE / "dna_tools.py"), "stats", "ghost.fasta"],
                       capture_output=True, text=True)
    assert r.returncode != 0 and "File not found" in (r.stderr + r.stdout)
    r = subprocess.run([sys.executable, str(HERE / "dna_tools.py"), "motif", str(demo)],
                       capture_output=True, text=True)
    assert r.returncode != 0  # missing motif argument fails cleanly
    print("PASS CLI end-to-end (stats, motif, missing file, missing motif)")


if __name__ == "__main__":
    tests = [v for k, v in sorted(globals().items()) if k.startswith("test_")]
    failed = 0
    for t in tests:
        try:
            t()
        except Exception:
            failed += 1
            import traceback
            print(f"FAIL {t.__name__}")
            traceback.print_exc()
    sys.exit(1 if failed else 0)
