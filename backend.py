from Bio import Blast
from Bio import SeqIO


def parser():
    for record in SeqIO.parse("name.fasta", "fasta"):
        sequenceID = record.id
        sequenceDescription = record.description
        sequence = str(record.seq)

