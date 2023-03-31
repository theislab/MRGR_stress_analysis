import argparse

parser = argparse.ArgumentParser(description='Diff Test for MRGR data')
parser.add_argument('-t', type=str,
                    help='Type of test',
                    default=100)

args = parser.parse_args()
print(args)

import os 
import numpy as np
import pandas as pd
import pickle
import scanpy as sc
import diffxpy.api as de
import os
os.chdir('./../')
print(os.getcwd())

adata = sc.read('../../data/processed/adata_reannotated_final.h5ad')
adata.obs['comb'] = adata.obs[['knockout', 'line', 'condition']].apply(lambda x: ' '.join(x), axis=1)
base_dir = '../../results/DE/de_stress/pickle/'

test = args.t

if test == 'WT':
    adata_loop = adata[adata.obs['line'] == 'Ctrl']
elif test == 'GR_Nex':
    adata_loop = adata[(adata.obs['line'] == 'Nex') & (adata.obs['knockout'] == 'GR')]
elif test == 'GR_Dlx':
    adata_loop = adata[(adata.obs['line'] == 'Dlx') & (adata.obs['knockout'] == 'GR')]
elif test == 'MR_Nex':
    adata_loop = adata[(adata.obs['line'] == 'Nex') & (adata.obs['knockout'] == 'MR')]
elif test == 'MR_Dlx':
    adata_loop = adata[(adata.obs['line'] == 'Dlx') & (adata.obs['knockout'] == 'MR')]
elif test == 'WT_GR':
    adata_loop = adata[(adata.obs['line'] == 'Ctrl') & (adata.obs['knockout'] == 'GR')]
elif test == 'WT_MR':
    adata_loop = adata[(adata.obs['line'] == 'Ctrl') & (adata.obs['knockout'] == 'MR')]

de_results = {'test': test}

coefs = de.utils.preview_coef_names(
sample_description=adata_loop.obs,
formula="~1+condition", 
)
           
for clust in adata_loop.obs['louvain_coarse'].cat.categories:
    adata_tmp = adata_loop[adata_loop.obs['louvain_coarse'] == clust,:].copy()
    adata_tmp.obs.size_factors = adata_tmp.obs.size_factors / np.mean(adata_tmp.obs.size_factors)
    print(f'In cluster {clust}:')
    sc.pp.filter_genes(adata_tmp, min_cells=np.floor(0.05*adata_tmp.shape[0]))
    print(f'Testing {adata_tmp.n_vars} genes...')
    print("")
    try:
        test_tmp = de.test.wald(
            data=adata_tmp.layers['counts'],
            formula_loc="~ 1 + condition + sample",
            size_factors='size_factors',
            factor_loc_totest='condition',
            sample_description=adata_tmp.obs,
            gene_names=adata_tmp.var_names,
            constraints_loc={'sample':'comb'},
            noise_model='nb',
            dtype="float64",
            )
        de_results[clust] = test_tmp
        de_results['coef_test'] = coefs[-1]
    except (np.linalg.LinAlgError, ValueError) as e:
        print('Model did not converge.')
        pass


with open(base_dir + 'de_results_' + test + '.pickle', 'wb') as handle:
    pickle.dump(de_results, handle, protocol=pickle.HIGHEST_PROTOCOL)
