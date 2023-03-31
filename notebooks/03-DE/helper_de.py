import scanpy as sc
import numpy as np
import pandas as pd
import pickle
from tqdm import tqdm
import scipy as sp
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
import matplotlib.colors as col
import seaborn as sns
import warnings

warnings.filterwarnings('ignore')

def de_postprocessing_new(
    de_results_listofdicts,
    adata,
    groupby='louvain_coarse',
    clean_sd=False,
    #clean_sd_2=False,
    save=None,
):
    de_post = []
    
    #Initial preprocessing
    for de_results_dict in de_results_listofdicts:
        test = de_results_dict['test']
        if test == 'WT':
            mask = ((adata.obs['line']=='Ctrl'))
        elif test == 'GR_Nex':
            mask = ((adata.obs['knockout'] == 'GR') & (adata.obs['line']=='Nex'))
        elif test == 'GR_Dlx':
            mask = ((adata.obs['knockout'] == 'GR') & (adata.obs['line']=='Dlx'))
        elif test == 'MR_Nex':
            mask = ((adata.obs['knockout'] == 'MR') & (adata.obs['line']=='Nex'))
        elif test == 'MR_Dlx':
            mask = ((adata.obs['knockout'] == 'MR') & (adata.obs['line']=='Dlx'))
        elif test == 'WT_GR':
            mask = ((adata.obs['line']=='Ctrl') & (adata.obs['knockout'] == 'GR'))
        elif test == 'WT_MR':
            mask = ((adata.obs['line']=='Ctrl') & (adata.obs['knockout'] == 'MR'))
        adata_test = adata[mask]
        
        print('Processing results for test {:s}'.format(de_results_dict['test']))
        for clust in de_results_dict:
            if clust not in ['test', 'coef_test']:
                tmp = de_results_dict[clust].summary()
                tmp = tmp.loc[tmp['ll'] != np.float64('-inf')]
                if clean_sd:
                    tmp = tmp.loc[tmp['coef_sd'] < 10]
                    tmp = tmp.loc[tmp['coef_sd'] > 0.01]
                #if clean_sd_2:
                #    tmp = tmp.loc[tmp['coef_sd'] < 10]
                #    tmp = tmp.loc[~(
                #        (tmp['coef_sd'] < 0.01)
                #        & (np.abs(tmp['log2fc']) > 10)
                #    )]
                #Sort values by q-value
                tmp.sort_values(by='qval', ascending=True, inplace=True)
                tmp = tmp.set_index('gene')
                tmp = tmp.assign(
                    clust = clust,
                    test = de_results_dict['test'],
                    mean_ctrl = np.nan,
                    mean_pert = np.nan
                )
                
                adata_clust = adata_test[adata_test.obs[groupby]==clust, adata_test.var_names.isin(tmp.index)].copy()
                adata_clust.var['mean_ctrl'] = np.mean(adata_clust.X[adata_clust.obs.condition=='No Stress'], axis=0)
                adata_clust.var['mean_pert'] = np.mean(adata_clust.X[adata_clust.obs.condition=='Stress'], axis=0)
                tmp['mean_pert'].loc[adata_clust.var_names] = adata_clust.var['mean_pert']
                tmp['mean_ctrl'].loc[adata_clust.var_names] = adata_clust.var['mean_ctrl']
                de_post.append(tmp)
            
    de_genes = pd.concat(de_post)
    de_genes.reset_index(inplace=True)
    if save:
        de_genes.to_pickle(save)
    return de_genes

def de_postprocessing_ctrl(
    de_results_listofdicts,
    adata,
    groupby='louvain_coarse',
    clean_sd=False,
    save=None,
):
    de_post = []
    
    #Initial preprocessing
    for de_results_dict in de_results_listofdicts:
        test = de_results_dict['test']
        if test == 'GR_Nex':
            mask = (
                (adata.obs['knockout'] == 'GR') 
                & (adata.obs['line'].isin(['Nex', 'Ctrl'])) 
                & (adata.obs['condition'] == 'No Stress')
            )
        elif test == 'GR_Dlx':
            mask = (
                (adata.obs['knockout'] == 'GR') 
                & (adata.obs['line'].isin(['Dlx', 'Ctrl'])) 
                & (adata.obs['condition'] == 'No Stress')
            )
        elif test == 'MR_Nex':
            mask = (
                (adata.obs['knockout'] == 'MR') 
                & (adata.obs['line'].isin(['Nex', 'Ctrl'])) 
                & (adata.obs['condition'] == 'No Stress')
            )
        elif test == 'MR_Dlx':
            mask = (
                (adata.obs['knockout'] == 'MR') 
                & (adata.obs['line'].isin(['Dlx', 'Ctrl'])) 
                & (adata.obs['condition'] == 'No Stress')
            )
        adata_test = adata[mask]
        
        print('Processing results for test {:s}'.format(de_results_dict['test']))
        for clust in de_results_dict:
            if clust not in ['test', 'coef_test']:
                tmp = de_results_dict[clust].summary()
                tmp = tmp.loc[tmp['ll'] != np.float64('-inf')]
                if clean_sd:
                    tmp = tmp.loc[np.abs(tmp['coef_sd']) > 0.01]
                    tmp = tmp.loc[np.abs(tmp['coef_sd']) < 10]
                #Sort values by q-value
                tmp.sort_values(by='qval', ascending=True, inplace=True)
                tmp = tmp.set_index('gene')
                tmp = tmp.assign(
                    clust = clust,
                    test = de_results_dict['test'],
                    mean_ctrl = np.nan,
                    mean_pert = np.nan
                )
                
                adata_clust = adata_test[adata_test.obs[groupby]==clust, adata_test.var_names.isin(tmp.index)].copy()
                adata_clust.var['mean_ctrl'] = np.mean(adata_clust.X[adata_clust.obs.line=='Ctrl'], axis=0)
                adata_clust.var['mean_pert'] = np.mean(adata_clust.X[adata_clust.obs.line!='Ctrl'], axis=0)
                tmp['mean_pert'].loc[adata_clust.var_names] = adata_clust.var['mean_pert']
                tmp['mean_ctrl'].loc[adata_clust.var_names] = adata_clust.var['mean_ctrl']
                de_post.append(tmp)
            
    de_genes = pd.concat(de_post)
    de_genes.reset_index(inplace=True)
    if save:
        de_genes.to_pickle(save)
    return de_genes


def tf_enrichment(adata, de_genes, mean_thr, log_thr, background_type='standard', direction=None):
    tf_genes_mr = pd.read_excel(
        'tf_targets.xlsx', 
        index_col=0, 
        sheet_name='MR-exclusive', 
        engine='openpyxl'
    )
    tf_genes_gr = pd.read_excel(
        'tf_targets.xlsx', 
        index_col=0, 
        sheet_name='GR-exclusive', 
        engine='openpyxl'
    )
    tf_genes_common = pd.read_excel(
        'tf_targets.xlsx', 
        index_col=0, 
        sheet_name='MR-GR-overlapping', 
        engine='openpyxl'
    )
    
    tf_genes_gr['tf'] = 'GR'
    tf_genes_mr['tf'] = 'MR'
    tf_genes_common_mr = tf_genes_common.copy()
    tf_genes_common_gr = tf_genes_common.copy()
    tf_genes_common_mr['tf'] = 'MR'
    tf_genes_common_gr['tf'] = 'GR'
    dataframes = [tf_genes_gr, tf_genes_common_gr, tf_genes_mr, tf_genes_common_mr]
    tf_genes = pd.concat(dataframes)
    tf_genes = tf_genes.rename(columns={"Gene symbol": 'target'})
    tf_genes = tf_genes.loc[tf_genes['target'].isin(adata.var_names)]
    
    nr3c1_set = np.array(tf_genes[tf_genes['tf']=='GR']['target'].values)
    nr3c2_set = np.array(tf_genes[tf_genes['tf']=='MR']['target'].values)
    
    de_genes = pd.read_pickle('./../../results/DE/de_stress_new/pickle/genes.pickle')
    de_genes_filter = de_genes[(de_genes['mean_pert']>mean_thr) | (de_genes['mean_ctrl']>mean_thr)]
    de_genes_filter = de_genes_filter[np.abs(de_genes_filter['log2fc'])>log_thr]

    de_genes_sign = de_genes[de_genes['qval']<0.05]
    de_genes_sign_mean = de_genes_sign[(de_genes_sign['mean_pert']>mean_thr) | (de_genes_sign['mean_ctrl']>mean_thr)]
    de_genes_sign_mean = de_genes_sign_mean[np.abs(de_genes_sign_mean['log2fc'])>log_thr]
    
    if direction=='up':
        de_genes_sign_mean = de_genes_sign_mean[de_genes_sign_mean['log2fc']>0]
    elif direction=='down':
        de_genes_sign_mean = de_genes_sign_mean[de_genes_sign_mean['log2fc']<0]

    
    listofdicts = []
    for test in ['WT', 'GR_Nex', 'GR_Dlx', 'MR_Nex', 'MR_Dlx']:
        print(test)
        for clust in adata.obs['louvain_coarse'].cat.categories:
            if background_type == 'standard':
                background = de_genes[
                            (de_genes['clust']==clust) 
                            & (de_genes['test']==test)
                        ].shape[0]
            elif background_type == 'filtered':
                background = de_genes_filter[
                            (de_genes_filter['clust']==clust) 
                            & (de_genes_filter['test']==test)
                        ].shape[0]

            tested_genes = de_genes_filter[
                        (de_genes_filter['clust']==clust) 
                        & (de_genes_filter['test']==test)
                    ]['gene'].values
            
            nr3c1_set_loop = set(nr3c1_set).intersection(set(tested_genes))
            nr3c2_set_loop = set(nr3c2_set).intersection(set(tested_genes))

            try:
                genes = de_genes_sign_mean[
                        (de_genes_sign_mean['clust']==clust) 
                        & (de_genes_sign_mean['test']==test)
                    ]

                sig_genes = genes['gene'].values

                inters = [
                len(set(sig_genes).intersection(set(nr3c1_set_loop))),
                len(set(sig_genes).intersection(set(nr3c2_set_loop)))
                ]
                enquiry = len(set(sig_genes))
                references = [
                    len(set(nr3c1_set_loop)),
                    len(set(nr3c2_set_loop))
                ]

                pval = hypergeom_test(
                        intersections=inters,
                        enquiry=enquiry,
                        references=references,
                        background=background
                )
                listofdicts.append({
                    'clust': clust,
                    'test': test,
                    'NDEGS': len(set(sig_genes)),
                    'GR overlap': inters[0],
                    'MR overlap': inters[1],
                    'GR pval': pval[0],
                    'MR pval': pval[1],
                    'GR ref': references[0],
                    'MR ref': references[1],
                    'background': background,
                    #'background_': background_,
                })
            except:
                pass

    tf_enr = pd.DataFrame(listofdicts)

    for i, row in tf_enr.T.iteritems():
        pvals = row[['GR pval', 'MR pval']].values
        adjusted_pvals = multipletests(pvals, alpha=0.05, method='fdr_bh')
        tf_enr.at[i, 'GR qval'] = adjusted_pvals[1][0]
        tf_enr.at[i, 'MR qval'] = adjusted_pvals[1][1]
    
        
    tf_enr['-log10 GR qval'] = -np.log10(tf_enr['GR qval'])
    tf_enr['-log10 MR qval'] = -np.log10(tf_enr['MR qval'])
    tf_enr = tf_enr[~tf_enr.test.isin(['WT_MR', 'WT_GR'])]
    tf_enr['clust'] = tf_enr['clust'].astype('category')
    tf_enr['clust'] = (
        tf_enr['clust']
        .cat
        .reorder_categories(
            adata.obs.louvain_coarse.cat.categories[
                adata.obs.louvain_coarse.cat.categories.isin(tf_enr.clust.unique())
            ]
        )
    )
    tf_enr['test'] = tf_enr['test'].astype('category')
    tf_enr['test'] = tf_enr['test'].cat.reorder_categories(
        ['WT', 'GR_Nex', 'GR_Dlx', 'MR_Nex', 'MR_Dlx']
    )
    return tf_enr
    
def plot_tf_enr(tf_enr, figsize=(10, 10)):
    
    cm_sign_gr = col.ListedColormap(['c'])
    fig, ax = plt.subplots(1, 2, figsize=figsize)
    d = tf_enr.pivot(index='clust', columns='test', values='-log10 GR qval')
    sns.heatmap(d, linewidths=.8, ax=ax[0], cmap='Blues', cbar=False)
    sns.heatmap(d[d>-np.log10(0.05)], linewidths=.8, ax=ax[0], cmap=cm_sign_gr, cbar=False)
    ax[0].set_title('GR enrichment')
    cm_sign_mr = col.ListedColormap(['r'])
    d = tf_enr.pivot(index='clust', columns='test', values='-log10 MR qval')
    sns.heatmap(d, linewidths=.8, ax=ax[1], cmap='Reds', cbar=False)
    sns.heatmap(d[d>-np.log10(0.05)], linewidths=.8, ax=ax[1], cmap=cm_sign_mr, cbar=False)
    ax[1].set_title('MR enrichment')
    
    fig.tight_layout()

def hypergeom_test(
        intersections: np.ndarray,
        enquiry: int,
        references: np.ndarray,
        background: int
) -> np.ndarray:
    """ Run a hypergeometric test (gene set enrichment).
    The scipy docs have a nice explanation of the hypergeometric ditribution and its parameters.
    This function wraps scipy.stats.hypergeom() and can compare multiple reference (enquiry)
    sets against one target set.
    :param intersections: np.ndarray
        Array with number of overlaps of reference sets with enquiry set.
    :param enquiry: np.ndarray
        Size of each enquiry set to be tested.
    :param references: int
        Array with size of reference sets.
    :param background: int
        Size of background set.
    """
    pvals = np.array([sp.stats.hypergeom(
        M=background,
        n=references[i],
        N=enquiry
    ).sf(x-1) for i, x in enumerate(intersections)])
    return pvals