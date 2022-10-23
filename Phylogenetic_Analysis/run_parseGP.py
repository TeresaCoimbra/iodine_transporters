import pandas as pd
from GPRecord import GPRecord
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO


def get_records(file_name):
    """Parse GenBank or GenPept file.
    Return records object"""
    records = SeqIO.parse(file_name,"gb")
    return records


#####

def parse_GP_recs(GenPept_records,output):
    seen_IDS = []  # List of viewed NCBI gene IDs
    seen_locus_tags = []  # list of viewed LocusTags
    acc = []  # list of Accessions
    acc_gene = { }
    my_records = []  # List of records to obtain a fasta file afterwards
    rows = []

    for r in GenPept_records:
        rec = GPRecord(r)
        gene = " "
        is_gene = True
        accession = rec.get_accession()
        source = rec.get_record_source()
        synonyms = rec.get_gene_syn()
        try:
            gene = rec.get_gene()
            is_gene = True
        except:
            try:
                gene = rec.get_LocusTag()
                is_gene = False
            except:
                gene = " "
                is_gene = False
            # except:
            # continue

        try:
            gene_ID = rec.get_NCBI_gene_id()
        except:
            gene_ID = ""

        other_ids = ""
        ### QUALIFIERS ###
        if rec.get_otherIDS() != "":
            other_ids = rec.get_otherIDS()

        try:
            defi = rec.get_record_description()
        except:
            continue

        if is_gene:
            if gene_ID not in seen_IDS:  # If the associated gene id wasn't seen yet, then create new entry

                new_row = { 'Accession': accession,
                            'Source': source,
                            'gene': gene,
                            'Definition': defi,
                            'Gene synonyms': synonyms,
                            'NCBI_gene_ID': gene_ID,
                            'DBS_ids': other_ids }
                seen_IDS.append(gene_ID)
                acc.append(accession)
                acc_gene[accession] = [source,gene]
                rows.append(new_row)

                rec_f = SeqRecord(rec.get_sequence(),id=accession,
                                  description=gene + " " + "[" + source + "]")  ### ADD SEQ TO FASTA
                my_records.append(rec_f)

        elif not is_gene:
            if gene not in seen_locus_tags:
                new_row = { 'Accession': accession,
                            'Source': source,
                            'gene': gene,
                            'Definition': defi,
                            'Gene synonyms': synonyms,
                            'NCBI_gene_ID': gene_ID,
                            'DBS_ids': other_ids }

                acc.append(accession)
                acc_gene[accession] = [source,gene]
                seen_locus_tags.append(gene)

                rows.append(new_row)
                rec_f = SeqRecord(rec.get_sequence(),id=accession,
                                  description=gene + " " + "[" + source + "]")
                my_records.append(rec_f)

    SeqIO.write(my_records,output,"fasta")
    df = pd.DataFrame.from_records(rows)
    print(df)
    return df,acc,acc_gene,my_records


########
def get_domains(GenPept_records,parsed_ids,accession_source_gene):
    """
    Checks for the existence of annotated conserved domains on a Genbank/GenPept file.

    Parameters
    ----------
    GenPept_records : Biopython object
        Created by parsing a GenBank or GenPept file
    parsed_ids      : List of non-redundant NCBI accessions
    accession_source_gene : dictionary {accession: [source, gene]}

    Returns
    -------
    domains
        list of dictionaries with domain regions
    """

    domains = []  # list with d

    for recx in GenPept_records:
        rec = GPRecord(recx)

        region = ''  # sequence region with domain: positions (x..y)
        r_name = ''  # name of the region
        note = ''  # associated notes
        CDD = ''  # CDD code

        accx = rec.get_accession()  # record accession

        if accx in parsed_ids:

            source = accession_source_gene[accx][0]
            gene = accession_source_gene[accx][1]

            for feat in range(len(recx.features)):  # uses recx instead of rec

                if recx.features[feat].type == "Region":
                    reg = recx.features[feat].location
                    region = str(reg).split("[")[1]
                    region = region.replace("]","").replace(":","..")

                    name = recx.features[feat].qualifiers["region_name"]
                    r_name = name[0]

                    try:
                        note = recx.features[feat].qualifiers["note"][0]
                    except:
                        note = " "

                    try:
                        CDD = recx.features[feat].qualifiers['db_xref'][0]
                    except:
                        CDD = " "

                entry = { 'Accession': accx,
                          'Source': source,
                          'gene': gene,
                          'Region': region,
                          'Region_name': r_name,
                          'CDD': CDD,
                          'other': note }

                if entry not in domains:
                    domains.append(entry)

    return domains


##################################################

def newFasta(df3,output):
    """Receives dataframe with the accessions that contain domain of interest
    Receives fasta file to filter - the result should be a fasta file with the selected accesions only"""
    new_recs = []
    records = SeqIO.parse(output,"fasta")
    print("*********************************************************")
    df4 = df3[["Accession","gene","Source"]].drop_duplicates()
    for seq_record in records:
        if seq_record.id in df4.Accession.values:
            new_recs.append(seq_record)

    output4 = output.replace(".faa","_new.faa")
    SeqIO.write(new_recs,output4,"fasta")
