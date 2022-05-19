from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd

# NAME OF GP FILE 
#fn = "danio_rerio_PDS"  # EXAMPLE


def qualifiers_to_dict(qualifier): 
    '''Takes a list of CDS qualifiers strings and returns dictionary with qualifiers as items'''
    sep = {}
    for i in qualifier:
        # i is a string 
        a = i.split(":") # split the string
        
        if a[0] == "HGNC" or a[0] =="MGI": # "HGNC appears as "HGNC:HGNC:code"
            sep[a[0]] = a[2]
        else: 
            sep[a[0]] = a[1]
    return sep
        

def parse_and_write(file_name):
    #df = create_dict()
    df = pd.DataFrame()
    
    records = SeqIO.parse(file_name, "gb")
    seen_IDS = []                                              # List of NCBI gene IDs
    seen_locus_tags = []
    acc = []                                                   # list of Accessions
    my_records = []                                            # List of records
    
    
    for record in records: 
        
        defi = record.description.split("[")                   # creates list regarding description
        source = defi[1].replace("]","")
        other_ids = ""
        
        try: 
            gene = record.features[-1].qualifiers["gene"][0]
        except:
            gene = record.features[-1].qualifiers["locus_tag"][0]
            
        ### QUALIFIERS ### 
        try:
            qualif = record.features[-1].qualifiers["db_xref"]     # qualifiers
            q = qualifiers_to_dict(qualif)
            gene_id = q["GeneID"]  
        
            
            for k, v in q.items():
                if (k != "GeneID") and (k != "CCDS"):
                    other_ids += str(k+':'+v+"; ")
        except:
            print("No gene ID for gene "+ gene)
        
        ### ADD ROW TO DF ###
        
        try:
            if gene_id not in seen_IDS:                            # If the associated gene id wasn't seen yet, then create new entry
        
                    new_row = {'Accession': record.id, 
                               'Source': source, 
                               'gene': gene,
                               'Definition': defi[0],
                               'NCBI_gene_ID': int(gene_id),
                               'DBS_ids': other_ids}
                    seen_IDS.append(gene_id)
                    acc.append(record.id)
                    seen_locus_tags.append(gene)

                    df = df.append(new_row, ignore_index=True)
                    rec = SeqRecord(record.seq, id = record.id, description = gene + " " + "["+ source + "]") ### ADD SEQ TO FASTA 
                    my_records.append(rec)
                    
        except:
            if gene not in seen_locus_tags:

                    new_row = {'Accession': record.id, 
                               'Source': source, 
                               'gene': gene,
                               'Definition': defi[0],
                               'NCBI_gene_ID': "",
                               'DBS_ids': other_ids}

                    acc.append(record.id)
                    seen_locus_tags.append(gene)
                    df = df.append(new_row, ignore_index=True)

                    rec = SeqRecord(record.seq, id = record.id, description = gene + " " + "["+ source + "]") ### ADD SEQ TO FASTA 
                    my_records.append(rec)

            

            
    SeqIO.write(my_records, output2, "fasta")
    
    return df, acc
  
  
  def parse_domains(file_name):
    '''Checks for the existence of conserved domains on a Genbank/GenPept file
    
    Returns
    -------
    Data frame with accession, source organism, name and domains. 
    '''
    
    records = SeqIO.parse(file_name, "gb")
    df, ids_conf = parse_and_write(file_name)
    
    #### CHECK DOMAINS ####
    df_domains = pd.DataFrame()
    for record in records: 
        region = ''
        r_name = ''
        note = ''
        CDD = ''

        if record.id in ids_conf:

            for feat in range(len(record.features)):

                if record.features[feat].type == "source":
                    s_org = record.features[feat].qualifiers["organism"][0]

                elif record.features[feat].type == "Region":
                    region = record.features[feat].location
                    region = str(region).split("[")[1]
                    region = region.replace("]","").replace(":","..")

                    r_name = record.features[feat].qualifiers["region_name"]
                    r_name = r_name[0]

                    note = record.features[feat].qualifiers["note"][0]

                    try: 
                        CDD = record.features[feat].qualifiers['db_xref'][0]
                    except:
                        continue

                try: 
                    n_gene = record.features[-1].qualifiers["gene"][0]

                except:
                    n_gene = record.features[-1].qualifiers["locus_tag"][0]

                new_row_dom = {'Accession': record.id,
                               'Source': s_org,
                               'gene': n_gene,
                               'Region': region,
                               'Region_name': r_name,
                               'CDD': CDD,
                               'other': note}

                if new_row_dom['Region'] != '':
                    df_domains = df_domains.append(new_row_dom, ignore_index=True)              
            
    return df_domains
  
#################################################################################################################
## PDS_files is a list of names of the downloaded GenPept files to process

for i in pds_files:
    fn = i 
    file = fn + ".gb"
    output1 = fn + ".xlsx"
    output2 = fn + ".faa"
    print("Output file names will be \n:", output1, "\n", output2)
    df1, ids_conf1 = parse_and_write(file)
    par_d = parse_domains(file)
    
    with pd.ExcelWriter(output1) as writer:  
        df1.to_excel(writer, sheet_name='info')
        par_d.to_excel(writer, sheet_name='domains')
    
    print("Number of sequences", len(df1.Source))
