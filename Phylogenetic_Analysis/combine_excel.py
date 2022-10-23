import os
import pandas as pd

### Combine excel files

data_folder = ""
output_file = "file.xlsx"

########

df_info = []
df_domains = []
for file in os.listdir(data_folder):
    if file.endswith('.xlsx'):
        print('Loading file {0}...'.format(file))
        df_info.append(pd.read_excel(os.path.join(data_folder, file), sheet_name="new_info"))
        df_domains.append(pd.read_excel(os.path.join(data_folder, file), sheet_name="domains"))

df_main = pd.concat(df_info, axis=0)
df_main2 = pd.concat(df_domains, axis=0)

with pd.ExcelWriter(output_file) as writer:
    df_main.to_excel(writer, sheet_name='info', index=False)
    df_main2.to_excel(writer, sheet_name='domains', index=False)
