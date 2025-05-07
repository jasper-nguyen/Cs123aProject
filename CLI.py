import argparse
from backend import load_sequence, viterbi, predict_starts

def main():
    parser = argparse.ArgumentParser(
        description="HMM‐based start‑site predictor for prokaryotic genes"
    )
    parser.add_argument("fasta", help="Input genome FASTA file")
    parser.add_argument("-o", "--output",
                        help="Write start positions to FILE (one per line)")
    args = parser.parse_args()

    seq = load_sequence(args.fasta)


    states = ["Intergenic", "Gene"]
    start_p = {"Intergenic": 0.95, "Gene": 0.05}
    trans_p = {
        "Intergenic": {"Intergenic": 0.999, "Gene": 0.001},
        "Gene":       {"Intergenic": 0.01,  "Gene": 0.99 }
    }
    emit_p = {
        "Intergenic": {"A":0.30, "C":0.20, "G":0.20, "T":0.30},
        "Gene":       {"A":0.20, "C":0.30, "G":0.30, "T":0.20},
    }

    path = viterbi(seq, states, start_p, trans_p, emit_p)
    starts = predict_starts(seq, path)

    if args.output:
        with open(args.output, "w") as out:
            for pos in starts:
                out.write(f"{pos}\n")
    else:
        for pos in starts:
            print(pos)

if __name__ == "__main__":
    main()