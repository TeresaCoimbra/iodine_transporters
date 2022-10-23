import os
import numpy as np
import pandas
from run_parseGP import *


def domains_to_df(domains_list):
    """
    domains_list : list of dictionary with annotated domains
    cond3 : third condition to add in order to select domains to add to dataframe

    Returns
    -------
        Dataframe with parsed domains
    """
    dm_rows = []

    for dom in domains_list:  # dom is a dictionary with domain information

        cond1 = dom['Region'] != ''
        cond2 = "Disordered" not in dom['Region_name']
        cond3 = "Anoctamin" in dom['Region_name']

        if (cond1) and (cond2) and cond3:
            dm_rows.append(dom)

    df_domains = pd.DataFrame.from_records(dm_rows)
    return df_domains


##################################
# Directory of gb files
path = r' '
os.chdir(path)

##################################

files = os.listdir(path)

pds_files = []  # list of GenPept files to be processed
for file in files:
    if ".gb" in file:
        pds_files.append(file)

##################################
for file_name in pds_files:
    global output2
    output1 = file_name.replace(".gb",".xlsx")  # create output name
    output2 = file_name.replace(".gb",".faa")

    print("Output file names will be: \n",output1,"\n",output2)
    print(file_name)
    # file = path + file_name

    ##########################################################################################
    ## Parse record info
    recs = get_records(file_name)  # Get records from a GenPept file
    df,acc,ASG,my_records = parse_GP_recs(recs,output2)  # Parse records from a GenPept file

    ##########################################################################################
    recs1 = get_records(file_name)
    dms = get_domains(recs1,acc,ASG)
    ###############################################
    domain_list = []
    df_dms = domains_to_df(dms)

    ######
    flag = True

    if "Taxus " in df.Source[0]:
        flag = False
    else:
        df_new = df.iloc[np.where(df.Accession.isin(df_dms.Accession))]
        print("*********************************************************")
        newFasta(df_new,output=output2)
    ######

    with pd.ExcelWriter(output1) as writer:
        df.to_excel(writer,sheet_name='info',index=False)
        df_dms.to_excel(writer,sheet_name='domains',index=False)
        if flag == True:
            df_new.to_excel(writer,sheet_name='new_info',index=False)
