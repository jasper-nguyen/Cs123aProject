from Bio import SeqIO
import math


#loads first sequence from the fasta and returns it as a string
def load_sequence(fasta_path: str) -> str:

    for rec in SeqIO.parse(fasta_path, "fasta"):
        return str(rec.seq).upper()
    raise ValueError(f"No sequences found in {fasta_path}")



#viterbi algorithm for hmm, gives sequence of parameters returning the most likely sequence
def viterbi(seq: str,
            states: list[str],
            start_p: dict[str, float],
            trans_p: dict[str, dict[str, float]],
            emit_p: dict[str, dict[str, float]]
           ) -> list[str]:

    V: list[dict[str, float]] = [{}]
    path: dict[str, list[str]] = {}

#calcualte probabilities for first base
    for s in states:
        e = emit_p[s].get(seq[0], 1e-8)
        V[0][s] = math.log(start_p[s]) + math.log(e)
        path[s] = [s]

#iterate over computing probabilities
    for t in range(1, len(seq)):
        V.append({})
        new_path: dict[str, list[str]] = {}
        for cur in states:
            best_prob, best_prev = max(         #finds best prevbious state that leads to current state
                (V[t-1][prev]
                 + math.log(trans_p[prev][cur])
                 + math.log(emit_p[cur].get(seq[t], 1e-8)),
                 prev)
                for prev in states
            )
            V[t][cur] = best_prob
            new_path[cur] = path[best_prev] + [cur]
        path = new_path

    final_state = max(states, key=lambda s: V[-1][s])
    return path[final_state]

#identifies translation start sites within gene region
def predict_starts(seq, state_path):
    starts = []
    for i in range(1, len(state_path) - 5):
        if state_path[i-1] == 'Intergenic' and state_path[i] == 'Gene':
            window = seq[i:i+6]
            if "ATG" in window:
                starts.append(i + window.index("ATG") + 1)
    return starts