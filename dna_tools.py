"""dna_tools.py - small DNA sequence utilities for FASTA files.

Reads a FASTA file and reports basic sequence statistics, reverse
complements, RNA transcriptions, or motif positions.

Usage:
  python3 dna_tools.py stats demo.fasta
  python3 dna_tools.py motif demo.fasta ATG
  python3 dna_tools.py transcribe demo.fasta
  python3 dna_tools.py revcomp demo.fasta

Sequences may contain A, C, G, T and N (unknown base), any case.
Pure Python standard library - no installs needed.
"""
import argparse
import sys
from pathlib import Path

VALID_BASES = set("ACGTN")


def read_fasta(path):
    """Parse a FASTA file into {name: sequence}. Raises ValueError on bad input."""
    records = {}
    name = None
    chunks = []
    with open(path) as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if name is not None:
                    records[name] = "".join(chunks)
                name = line[1:].strip() or "unnamed"
                if name in records:
                    raise ValueError(f"Duplicate sequence name: {name}")
                chunks = []
            else:
                if name is None:
                    raise ValueError("Found sequence data before the first '>' header - not a FASTA file.")
                chunks.append(line.upper())
    if name is not None:
        records[name] = "".join(chunks)
    if not records:
        raise ValueError(f"No sequences found in {path}.")
    for n, seq in records.items():
        if not seq:
            raise ValueError(f"Sequence '{n}' is empty.")
        bad = set(seq) - VALID_BASES
        if bad:
            raise ValueError(f"Sequence '{n}' contains invalid base(s): {', '.join(sorted(bad))}. Only A, C, G, T, N allowed.")
    return records


def gc_content(seq):
    """GC fraction (0-100) over known bases only; Ns are excluded."""
    known = [b for b in seq if b != "N"]
    if not known:
        raise ValueError("Sequence has no known bases (all N).")
    return round(100.0 * sum(1 for b in known if b in "GC") / len(known), 1)


def reverse_complement(seq):
    comp = {"A": "T", "T": "A", "C": "G", "G": "C", "N": "N"}
    return "".join(comp[b] for b in reversed(seq))


def transcribe(seq):
    """DNA to RNA (T -> U). N stays N."""
    return seq.replace("T", "U")


def find_motif(seq, motif):
    """1-based start positions of every occurrence, including overlaps."""
    motif = motif.upper()
    if not motif:
        raise ValueError("Motif must not be empty.")
    bad = set(motif) - VALID_BASES
    if bad:
        raise ValueError(f"Motif contains invalid base(s): {', '.join(sorted(bad))}.")
    hits = []
    start = seq.find(motif)
    while start != -1:
        hits.append(start + 1)
        start = seq.find(motif, start + 1)
    return hits


def main(argv=None):
    p = argparse.ArgumentParser(description="DNA sequence utilities for FASTA files.")
    p.add_argument("command", choices=["stats", "motif", "transcribe", "revcomp"])
    p.add_argument("fasta", help="input FASTA file")
    p.add_argument("motif", nargs="?", help="motif to search for (motif command only)")
    args = p.parse_args(argv)

    src = Path(args.fasta)
    if not src.exists():
        sys.exit(f"File not found: {src}")
    try:
        records = read_fasta(src)
    except ValueError as e:
        sys.exit(f"Invalid input: {e}")

    if args.command == "motif" and not args.motif:
        sys.exit("The motif command needs a motif, e.g.: python3 dna_tools.py motif demo.fasta ATG")
    if args.command != "motif" and args.motif:
        sys.exit(f"The {args.command} command does not take a motif argument.")

    for name, seq in records.items():
        if args.command == "stats":
            print(f"{name}: length {len(seq)}, GC {gc_content(seq)}%")
        elif args.command == "revcomp":
            print(f">{name} reverse-complement\n{reverse_complement(seq)}")
        elif args.command == "transcribe":
            print(f">{name} transcribed (RNA)\n{transcribe(seq)}")
        elif args.command == "motif":
            hits = find_motif(seq, args.motif)
            if hits:
                print(f"{name}: motif {args.motif.upper()} found at position(s) {', '.join(map(str, hits))}")
            else:
                print(f"{name}: motif {args.motif.upper()} not found")
    return 0


if __name__ == "__main__":
    sys.exit(main())
