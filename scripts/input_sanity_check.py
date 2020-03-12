from Bio import SeqIO
from Bio import SearchIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from scripts.load_mappings import load_sequences

def check_sequences(filename, out_dir):
    seqs = load_sequences(filename, "list")
    unique_seqs = []
    for i in range(len(seqs)):
        add = True
        for j in range(i + 1, len(seqs)):
            if seqs[i].seq == seqs[j].seq or seqs[i].id == seqs[j].id:
                add = False
                break
        if add:
            unique_seqs.append(seqs[i])
    if len(unique_seqs) < len(seqs):
        save_fasta(out_dir + "/unique_sequences", unique_seqs)
        return False, out_dir + "/unique_sequences.fasta"
    else:
        return True, filename

def save_fasta(filename, orfs):
    with open(filename + ".fasta", "w") as output_handle:
        SeqIO.write(orfs, output_handle, "fasta")