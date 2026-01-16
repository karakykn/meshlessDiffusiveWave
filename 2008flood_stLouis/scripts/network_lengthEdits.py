import pandas as pd
import numpy as np

for nn in range(5):
    if nn == 0:
        df = pd.read_csv(f'../data/from_saverton_mississippi.csv')
    elif nn == 1:
        df = pd.read_csv(f'../data/from_hermann_missouri.csv')
    elif nn == 2:
        df = pd.read_csv(f'../data/from_valleycity_illinois.csv')
    elif nn == 3:
        df = pd.read_csv(f'../data/from_grafton_mississippi.csv')
    else:
        df = pd.read_csv(f'../data/from_madison_mississippi.csv')

    nodeN = df.shape[0]
    remover = []
    for i in range(nodeN-1):
        if df['Length'].iloc[i] < 4000:
            remover.append(i)
            df['Length'].iloc[i+1] += df['Length'].iloc[i]
    # df = df.drop(remover)
    if df['Length'].iloc[-1] < 4000:
        df['Length'].iloc[i - 1] += df['Length'].iloc[i]
        remover.append(nodeN-1)
    df = df.drop(remover)

    nodeN = df.shape[0]

    for ret in range(5):
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

    if nn == 0:
        df.to_csv(f'../data/from_saverton_mississippi_edit.csv', index = False)
    elif nn == 1:
        df.to_csv(f'../data/from_hermann_missouri_edit.csv', index = False)
    elif nn == 2:
        df.to_csv(f'../data/from_valleycity_illinois_edit.csv', index = False)
    elif nn == 3:
        df.to_csv(f'../data/from_grafton_mississippi_edit.csv', index = False)
    else:
        df.to_csv(f'../data/from_madison_mississippi_edit.csv', index = False)