# From a fasta file obtain only entries for the selected species
import re
from Bio import SeqIO


########################################################
def removing_ids_fasta(fasta_f,out_file,selected_sps):
    recs = []
    for seq_record in SeqIO.parse(fasta_f,"fasta"):
        pattern = "\[(.*?)\]"
        organism = re.search(pattern,seq_record.description).group(1)
        if organism in selected_sps:
            recs.append(seq_record)
        else:
            continue

    SeqIO.write(recs,out_file,"fasta")

    print("Processing complete")


########################################################
# Sps are the species to select
sps = ["Homo sapiens","Mus musculus","Danio rerio","Ciona intestinalis","Octopus bimaculoides",
       'Hydra vulgaris',"Amphimedon queenslandica","Monosiga brevicollis MX1","Dictyostelium discoideum AX4",
       "Saccharomyces cerevisiae S288C","Candida albicans SC5314","Schizosaccharomyces pombe",
       "Arabidopsis thaliana","Oryza sativa Japonica group","Amborella trichopoda",
       "Selaginella moellendorffii","Chlamydomonas reinhardtii","Micromonas commoda"]

fasta_f = "file.fasta"
output_file = "file2.fasta"

removing_ids_fasta(fasta_f,output_file,sps)
