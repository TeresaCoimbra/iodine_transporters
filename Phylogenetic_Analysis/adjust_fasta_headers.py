### Adjust header of a fasta file from ">Accession Gene [organism]" to ">(ORGANISM abbreviation) Gene | Acession"
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
import re

fasta_file = r'file.fasta'
output_file = fasta_file.replace(".fasta","_header.fasta")



def change_header(fasta_f, output):
        abb = { 'Amphimedon queenslandica': 'Aq','Arabidopsis thaliana': 'At',
                'Aspergillus niger CBS 513.88': 'Na',
                'Amborella trichopoda': 'Atr','Candida albicans SC5314': 'Ca',
                'Chondrus crispus': 'Cc','Caenorhabditis elegans': 'Ce',
                'Ciona intestinalis': 'Ci','Cryptococcus neoformans var. neoformans JEC21': 'Cn',
                'Chlamydomonas reinhardtii': 'Cr','Dictyostelium discoideum AX4': 'Dd',
                'Drosophila melanogaster': 'Dm','Danio rerio': 'Dr',
                'Ectocarpus siliculosus': 'Es','Homo sapiens': 'Hs',
                'Hydra vulgaris': 'Hv','Monosiga brevicollis MX1': 'Mb',
                'Micromonas commoda': 'Mc','Mus musculus': 'Mm','Neurospora crassa OR74A': 'Nc',
                'Octopus bimaculoides': 'Ob','Oryza sativa Japonica Group': 'Os',
                'Physcomitrella patens ': 'Pp','Physcomitrium patens': 'Pp',
                'Saccharomyces cerevisiae S288C': 'Sc','Selaginella moellendorffii': 'Sm',
                'Schizosaccharomyces pombe': 'Sp','Taxus chinensis': 'Tc',
                'Ustilago maydis 521': 'Um','Yarrowia lipolytica CLIB122': 'Yl' }

        pattern = "\[(.*?)\]"
        pattern2 = "\s(.*?)\s"
        all_records = []
        for seq_record in SeqIO.parse(fasta_f,"fasta"):
                species = re.search(pattern,seq_record.description).group(1)
                if species in abb.keys():
                        source = abb[species]
                accession = seq_record.id
                gene_name = re.search(pattern2,seq_record.description).group(1)
                rx = SeqRecord(seq_record.seq, id = "("+source+") "+ gene_name, description = "| " + accession)
                all_records.append(rx)

        SeqIO.write(all_records, output, "fasta")


change_header(fasta_file,output_file)
