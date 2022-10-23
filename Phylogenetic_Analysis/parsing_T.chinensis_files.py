# Special parsing of GB files from Taxus chinensis (incomplete annotations)
from run_parseGP import *

file_name = " " # Taxus GB file
output1 = file_name.replace(".gb", ".xlsx")
output2 = file_name.replace(".gb", ".faa")

recs = get_records(file_name)  # Get records from a GenPept file
df, acc, ASG, my_records = parse_GP_recs(recs, output2)  # Parse records from a GenPept file

with pd.ExcelWriter(output1) as writer:
    df.to_excel(writer,sheet_name='info',index=False)
