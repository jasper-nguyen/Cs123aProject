#!/usr/bin/env python3
import argparse
from backend import load_sequence, viterbi, predict_starts

def main():
    print("HMM-based start-site predictor\n")


    parser = argparse.ArgumentParser(
        description="HMM-based start-site predictor for prokaryotic genes"
    )
    parser.add_argument("fasta", help="Input genome FASTA file")
    parser.add_argument("-o", "--output",
                        help="Write predicted start positions to FILE (one per line)")
    args = parser.parse_args()


    seq = load_sequence(args.fasta)
    print(f"Loaded sequence: {len(seq):,} bp from '{args.fasta}'")

    all_atgs = [i+1 for i in range(len(seq) - 2) if seq[i:i+3] == 'ATG']
    print(f"ATG codons detected: {len(all_atgs)} (first 10 at positions: {all_atgs[:10]})")


    states = ["Intergenic", "Gene"]
    start_p = {"Intergenic": 0.85, "Gene": 0.15}
    trans_p = {
        "Intergenic": {"Intergenic": 0.90, "Gene": 0.10},
        "Gene": {"Intergenic": 0.10, "Gene": 0.90}
    }
    emit_p = {
        "Intergenic": {"A": 0.25, "C": 0.25, "G": 0.25, "T": 0.25},
        "Gene": {"A": 0.60, "C": 0.05, "G": 0.30, "T": 0.05}
    }


    path = viterbi(seq, states, start_p, trans_p, emit_p)
    print(f"HMM path computed across {len(path):,} states")
    print(f"Total positions in 'Gene' state: {path.count('Gene'):,}")


    starts = predict_starts(seq, path)
    print(f"Predicted start sites: {len(starts)}")
    if starts:
        print("Start site positions:", starts[:10])
    else:
        print("No start sites predicted.")


    if args.output:
        with open(args.output, "w") as out:
            for pos in starts:
                out.write(f"{pos}\n")

if __name__ == "__main__":
    main()
