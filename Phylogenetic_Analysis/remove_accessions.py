
import pandas as pd

# remove accessions from a dataframe (or excel file)


def read_txt_ids(file_name):
    '''Reads txt files with ids to be removed, considering that IDs have been written one per line.

    Returns
    -------
        List with ids without space or paragraph marks'''

    f = open(file_name, "r")
    content_list = f.readlines()
    new_list = []
    for i in content_list:
        new_list.append(i.replace("\n", "").replace(" ", ""))

    return new_list

def remove_acx(df,acx):
    '''Remove rows of dataframe. Dataframe must have a column with accessions.
    Accessions to remove are listed in a txt file (ids)'''

    reslt = df[~df["Accession"].isin(acx)]
    return reslt


### run
# input definition
data_file = r'data.xlsx'
txt_file = r'data.txt'

# output definition
output_file = r'data.xlsx'

###
df1 = pd.read_excel(data_file, sheet_name="info")
df2 = pd.read_excel(data_file, sheet_name="domains")

ids = read_txt_ids(txt_file)   # must be a list
df1_1 = remove_acx(df1, ids)
df2_2 = remove_acx(df2, ids)


flag = input("Save to excel? Y/N: ")
if flag == 'Y':
    with pd.ExcelWriter(output_file) as writer:
        df1_1.to_excel(writer, sheet_name='info', index=False)
        df2_2.to_excel(writer, sheet_name='domains', index=False)
else:
    print("Dataset not saved.")
    print("Species", df1_1.Source.unique())
    print("Total of", len(df1_1.Source.unique()), "species")
