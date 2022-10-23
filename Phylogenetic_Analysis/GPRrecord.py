from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pandas as pd
import re


def qualifiers_to_dict(qualifier):
    '''Takes a list of CDS qualifiers strings and returns dictionary with qualifiers as items'''
    sep = {}
    for i in qualifier:
        # i is a string
        a = i.split(":")  # split the string

        if a[0] == "HGNC" or a[0] == "MGI":  # "HGNC appears as "HGNC:HGNC:code"
            sep[a[0]] = a[2]
        else:
            sep[a[0]] = a[1]
    return sep


class GPRecord:

    def __init__(self, record):
        self.record = record
        self.gene_id = ""
        self.gene = ""
        self.locus_tag = ""

    def get_NCBI_gene_id(self):
        qualif = self.record.features[-1].qualifiers["db_xref"]  # qualifiers
        q = qualifiers_to_dict(qualif)
        self.gene_id = q["GeneID"]
        return self.gene_id

    def get_record_description(self):
        '''From a record description obtain description of sequence'''
        defi = self.record.description.split("[")  # creates list regarding description
        self.description = defi[0]
        return self.description

    def get_record_source(self):
        '''Get source organism of record in GenPept file'''
        defi = self.record.description
        pattern = "\[(.*?)\]"
        try:
            self.source = re.search(pattern, defi).group(1)
        except:
            self.source = "Undefined"
        return self.source

    def get_accession(self):
        self.accession = self.record.id
        return self.accession

    def get_gene(self):
        self.gene = self.record.features[-1].qualifiers["gene"][0]
        return self.gene

    def get_LocusTag(self):
        self.locus_tag = self.record.features[-1].qualifiers["locus_tag"][0]
        return self.locus_tag

    def get_otherIDS(self):
        '''Returns string with other IDS from external databases, if present'''
        self.other_ids = ''
        try:
            qualif = self.record.features[-1].qualifiers["db_xref"]  # qualifiers
            q = qualifiers_to_dict(qualif)
            for k, v in q.items():
                if (k != "GeneID") and (k != "CCDS"):
                    self.other_ids += str(k + ':' + v + "; ")
        except:
            print(" ")
        return self.other_ids


    def get_gene_syn(self):
        '''Gets the gene synonyms, if any are anotated in the GB file'''
        self.synonyms = ''
        try:
            self.synonyms = self.record.features[-1].qualifiers["gene_synonym"]  # qualifiers
        except:
            self.synonyms = ''

        return self.synonyms

    def get_sequence(self):
        return self.record.seq

