# DNA Tools

Small command-line utilities for DNA sequences in FASTA files: sequence
statistics (length, GC content), reverse complement, RNA transcription, and
motif search with 1-based positions (overlapping matches included).

Pure Python standard library - no installs, no accounts, runs locally.

## Why it exists

A portfolio project joining biology and code: correct FASTA parsing with
strict validation (A, C, G, T, N only), and sequence operations implemented
from scratch rather than imported from a bioinformatics library. `demo.fasta`
contains three short made-up sequences written for this project, so it works
out of the box.

## Usage

```bash
python3 dna_tools.py stats demo.fasta          # length and GC% per sequence
python3 dna_tools.py motif demo.fasta ATG      # 1-based positions of every ATG
python3 dna_tools.py transcribe demo.fasta     # DNA -> RNA
python3 dna_tools.py revcomp demo.fasta        # reverse complement
```

Input is validated: invalid bases, missing headers, empty files and missing
arguments all fail with a clear message instead of a traceback.

## Tests

```bash
python3 test_dna_tools.py
```

Six tests covering the GC% maths (including N exclusion), reverse complement,
transcription, motif positions (including overlapping matches), FASTA parsing
errors, and a full command-line run.

## Scope and limitations

This is a small teaching-scale tool, not a replacement for real
bioinformatics software. It handles only DNA (no protein sequences, no
quality scores), reads whole files into memory (fine for small FASTA files,
not genomes), and does exact motif matching only - no mismatches, gaps or
regular-expression motifs. N bases are supported but excluded from the GC%
denominator.

*Author: Varun Suresh*
