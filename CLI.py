#!/usr/bin/env python3
import argparse
from backend import load_sequence, viterbi, predict_starts

def main():
    # ─── Debug #1: show we’re in main() ────────────────────────
    print(">>> CLI starting…")

    parser = argparse.ArgumentParser(
        description="HMM‐based start‑site predictor for prokaryotic genes"
    )
    parser.add_argument("fasta", help="Input genome FASTA file")
    parser.add_argument("-o", "--output",
                        help="Write start positions to FILE (one per line)")
    args = parser.parse_args()

    # ─── Debug #2: confirm sequence length ─────────────────────
    seq = load_sequence(args.fasta)
    print(f">>> Loaded {len(seq)} bp from {args.fasta}")

    # Optional extra check: do we even have any ATGs?
    all_atgs = [i+1 for i in range(len(seq)-2) if seq[i:i+3] == 'ATG']
    print(f">>> Found {len(all_atgs)} total ATG codons (first 10): {all_atgs[:10]}")

    # --- HMM parameters (defaults) ---
    states = ["Intergenic", "Gene"]
    start_p = {"Intergenic": 0.85, "Gene": 0.15}
    trans_p = {
        "Intergenic": {"Intergenic": 0.90, "Gene": 0.10},
        "Gene": {"Intergenic": 0.10, "Gene": 0.90}
    }

    # emission odds that favor ATG in the Gene state:
    emit_p = {
        "Intergenic": {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
        "Gene": {"A": 0.6, "C": 0.05, "G": 0.3, "T": 0.05}
    }

    # ─── Debug #3: run the HMM ─────────────────────────────────
    path = viterbi(seq, states, start_p, trans_p, emit_p)
    print(f">>> HMM path computed ({len(path)} states)")
    print(f">>> Gene state count: {path.count('Gene')}")
    # ─── Debug #4: scan for start‐codon transitions ─────────────
    starts = predict_starts(seq, path)
    print(f">>> predict_starts found {len(starts)} candidate sites")

    # Output
    if args.output:
        with open(args.output, "w") as out:
            for pos in starts:
                out.write(f"{pos}\n")
    else:
        for pos in starts:
            print(pos)

    # ─── Debug #5: print first 10 starts for sanity ────────────
    print(">>> first 10 starts:", starts[:10])


if __name__ == "__main__":
    main()
