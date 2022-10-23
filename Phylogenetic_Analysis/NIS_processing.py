import os
import pandas as pd
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
    print(domains_list)
    print('*****')

    for dom in domains_list:  # dom is a dictionary with domain information
        # print(dom)
        cond1 = dom['Region'] != ''
        cond2 = "Disordered" not in dom['Region_name']

        mfs_dms = ['CDD:271393','CDD:271402','CDD:212089','CDD:271394','CDD:271393',
                   'CDD:271402','CDD:271393','CDD:271383','CDD:271382','CDD:271379',
                   'CDD:271381','CDD:212056','CDD:271380','CDD:212058','CDD:271368',
                   'CDD:382020','CDD:271370']
        cond3 = dom['CDD'] in mfs_dms

        if cond1 and cond2 and cond3:
            dm_rows.append(dom)

    df_domains = pd.DataFrame.from_records(dm_rows)
    return df_domains


##############################################################################
# Directory of gb files

path = r''
print("Current path is: ",path)
a = input("Use another path to run? Y/N: ")
if a == "Y":
    path = input("Insert alternative path to folder containing gb files: ")
os.chdir(path)
# DIRECT = input("Directory of folder with GenPept files to process:")

##################################

files = os.listdir(path)

pds_files = []  # list of GenPept files to be processed
for file in files:
    if ".gb" in file:
        pds_files.append(file)

##################################
for file_name in pds_files:
    output1 = file_name.replace(".gb",".xlsx")  # create output name
    output2 = file_name.replace(".gb",".faa")
    output3 = file_name.replace(".gb","_new.faa")

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

    if "Taxus" in df.Source[0]:
        flag = False
    else:
        try:
            df_new = df.iloc[np.where(df.Accession.isin(df_dms.Accession))]
            newFasta(df_new,output=output2)
        except:
            print("No filter possible for ",df["Source"][0])

    ######

    with pd.ExcelWriter(output1) as writer:
        df.to_excel(writer,sheet_name='info',index=False)
        df_dms.to_excel(writer,sheet_name='domains',index=False)
        df_new.to_excel(writer,sheet_name='new_info',index=False)
