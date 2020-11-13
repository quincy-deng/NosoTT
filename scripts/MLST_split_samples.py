import os
import pandas as pd
import sys


def quit_with_error(message):
    """Displays the given message and ends the program's execution."""
    print('Error:', message, file=sys.stderr)
    sys.exit(1)


def get_samples(names, casedate, items):
    if os.path.exists('data/mlst_database/STs.txt'):
        df = pd.read_csv('data/mlst_database/STs.txt', sep='\t',
                         names=['name', 'species', 'ST'] + ["gene" + str(i) for i in range(7)])
        fastas = [i for i in names if i not in list(df['name'])]
    if len(fastas) > 0:
        os.system('mlst {} >> data/mlst_database/STs.txt'.format(' '.join(fastas)))
    df = pd.read_csv('data/mlst_database/STs.txt', sep='\t',
                     names=['name', 'species', 'ST'] + ["gene" + str(i) for i in range(7)])
    df_temp = pd.DataFrame()
    df_temp['item'] = items
    df_temp['name'] = names
    df = pd.merge(df_temp, df, how='left', on='name')
    df = df[~df['ST'].isin(['-'])]
    df['spe_st'] = df['species'].str.cat(df['ST'].apply(str), sep="_ST")
    df_merge = df.groupby(by='spe_st').apply(lambda x: '_'.join(x['item'].apply(str)))
    df_merge2 = df.groupby(by='spe_st').apply(lambda x: '_'.join(x['name'].apply(str)))
    index = df_merge.index
    name = [(k1, k2) for k1, k2 in zip(df_merge, df_merge2) if len(k1.split('_')) > 1]
    x = {'.'.join([i, j]): k.split('_') for i, (j, k) in zip(index, name)}
    return x


def get_sample_config(fastas, samples):
    df = pd.read_csv('data/casedate/temp_df.csv')
    x = [i for i in samples if i not in map(str, df['item'])]
    if len(x) != 0:
        quit_with_error("sample {} has no clinical data".format(",".join(x)))

    y = [os.path.basename(i) for i in fastas if not os.path.exists(i)]
    if len(y) != 0:
        quit_with_error("sample {} has no fasta file".format(",".join(y)))

    samples_items = get_samples(fastas, 'data/casedate/casedatabase.csv', list(df['item']))
    return samples_items
