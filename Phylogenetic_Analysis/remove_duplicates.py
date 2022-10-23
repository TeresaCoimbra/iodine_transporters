### Removing duplicate sequences from a fasta file

from Bio import SeqIO
######################################################
file_in = r'exp1.fasta'

######################################################


def remove_duplicates(fasta_file):

    assert ".fasta" in fasta_file, "Input file is not a fasta file"
    seen = set()
    records = []

    for record in SeqIO.parse(fasta_file, "fasta"):
        if record.seq not in seen:
            seen.add(record.seq)
            records.append(record)

    file_out = fasta_file.replace(".fasta", "_no_dup.fasta")
    #writing to a new fasta file
    SeqIO.write(records, file_out, "fasta")
    print("done")

remove_duplicates(file_in)
