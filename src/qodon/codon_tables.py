import python_codon_tables as pct
import pandas as pd
import numpy as np

def convert_to_p_list(a):
    '''
    Helper function

    '''
    j = 0
    l = []
    for i in a:
        j += i
        l.append(j)
    return l


def construct_codon_table(species="e_coli_316407"):
    '''
    Build reference table containing:

        amino acid:codon mappings
        codon frequencies

    This data is referenced by both the GA and the BQM

    '''
    # Load codon data
    codons_tables = pct.get_all_available_codons_tables()
    table = pct.get_codons_table(species)
    df = pd.DataFrame([(a, c, s) for a, v in table.items()
                       for c, s in v.items() if a != '*'],
                      columns=['aa', 'codon', 'score'])

    # Transform data into useful format
    df['tup'] = df.apply(lambda x: (x['codon'], x['score']), axis=1)
    by_aa = df.groupby('aa')
    ms_by_aa = by_aa['tup'].apply(list).apply(
        lambda x: max(x, key=lambda l: l[1]))
    df['log_score'] = df.apply(
        lambda x: abs(np.log(x['score'] / ms_by_aa.loc[x['aa']][1])), axis=1)

    # Merge lists of data into dataframe
    code_map_2 = pd.DataFrame(by_aa['score'].apply(list))
    code_map_2 = code_map_2.merge(pd.DataFrame(by_aa['codon'].apply(list)),
                                  left_index=True,
                                  right_index=True)
    code_map_2 = code_map_2.merge(pd.DataFrame(by_aa['log_score'].apply(list)),
                                  left_index=True,
                                  right_index=True)
    code_map_2.rename(columns={
        'score': 'scores',
        'codon': 'codons',
        'log_score': 'log_scores'
    },
                      inplace=True)

    # Convert dataframe to dict for quick lookups
    code_map_2['probs'] = code_map_2['scores'].apply(convert_to_p_list)
    code_map = code_map_2.to_dict('index')
    codon_scores = dict([
        item for sublist in
        [list(zip(_['codons'], _['log_scores'])) for _ in code_map.values()]
        for item in sublist
    ])
    codon_scores = {k: abs(v) for k, v in codon_scores.items()}

    codon_table = {k: v['codons'] for k, v in code_map.items()}

    return codon_table, codon_scores, code_map


codon_table, codon_scores, code_map = construct_codon_table()
