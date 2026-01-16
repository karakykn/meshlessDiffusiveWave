import pandas as pd
import numpy as np

df = pd.read_csv(f'../data/from_tarbert.csv')

nodeN = df.shape[0]
remover = []
for i in range(nodeN-1):
    if df['Length'].iloc[i] < 1000:
        remover.append(i)
        df['Length'].iloc[i+1] += df['Length'].iloc[i]
df = df.drop(remover)

nodeN = df.shape[0]

for _ in range(5):
    i = 0
    while i < nodeN:
        if df['Length'].iloc[i] > 9999:
            new_row = df.iloc[i]
            new_row['Length'] = df['Length'].iloc[i] / 2
            df['Length'].iloc[i] = df['Length'].iloc[i] / 2
            df = pd.concat(
                [df.iloc[:i], pd.DataFrame([new_row]), df.iloc[i:]]
            ).reset_index(drop=True)
            nodeN += 1
        i += 1

    # mask = df['Length'] > 9999  # condition
    # rows_to_dup = df[mask]
    # df = pd.concat(
    #     [df.iloc[:index], pd.DataFrame([new_row]), df.iloc[index:]]
    # ).reset_index(drop=True)

df.to_csv(f'../data/fromTarbert_edit.csv')