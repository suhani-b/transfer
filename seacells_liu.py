# FILE: seacells_script.py
# AUTHOR: Suhani Balachandran, Sam Rose
# DATE: June 2024
# DESCRIPTION: Runs SEACells on a specified sample in a dataset
# PARAMS: path_to_anndata, study_name, sample_name, path_to_results_dir, {optional} lineage_key=""
# RESULT: sample-specific directory within the specified results_dir that contains
#     a folder of plots, a subsetted and modified anndata object, and summary statistics 
# EXAMPLE RUN: python3 seacells_script.py /full_anndata_path Salcher He_Fan_2021_LUAD1 seacells_results {lineage_key}

import numpy as np
import pandas as pd
import scanpy as sc
import SEACells
import matplotlib as mpl
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import json
from SEACells.core import summarize_by_SEACell
import palantir
from scipy import stats
import os

# ----------------------------------------------------------------------
# HELPER FUNCTIONS
# plotting
def plot_umap(tmp_ad, cell_type_key):
    fig = sc.pl.umap(tmp_ad, color=cell_type_key, return_fig=True, size=10)
    h, l = fig.axes[0].get_legend_handles_labels()
    color_values = [(0.8392156862745098, 0.15294117647058825, 0.1568627450980392), (0.7372549019607844, 0.7411764705882353, 0.13333333333333333), (1.0, 0.4980392156862745, 0.054901960784313725), (0.5803921568627451, 0.403921568627451, 0.7411764705882353), (0.12156862745098039, 0.4666666666666667, 0.7058823529411765), (0.17254901960784313, 0.6274509803921569, 0.17254901960784313), (0.5490196078431373, 0.33725490196078434, 0.29411764705882354), (0.7803921568627451, 0.7803921568627451, 0.7803921568627451), (0.8901960784313725, 0.4666666666666667, 0.7607843137254902), (0.596078431372549, 0.8745098039215686, 0.5411764705882353), (0.6823529411764706, 0.7803921568627451, 0.9098039215686274), (1.0, 0.596078431372549, 0.5882352941176471), (0.4980392156862745, 0.4980392156862745, 0.4980392156862745),  (0.8588235294117647, 0.8588235294117647, 0.5529411764705883), (0.7725490196078432, 0.6901960784313725, 0.8352941176470589), (0.7686274509803922, 0.611764705882353, 0.5803921568627451), 'cornsilk', (0.9686274509803922, 0.7137254901960784, 0.8235294117647058)]
    color_dict = dict([(l[x], color_values[x]) for x in range(len(l))])

    df = pd.DataFrame(tmp_ad.obsm['X_umap']).set_index(tmp_ad.obs_names).join(tmp_ad.obs['SEACell']).groupby('SEACell').mean()
    df = df.join(tmp_ad.obs[['SEACell',cell_type_key]].groupby('SEACell').agg(lambda x:x.value_counts().index[0]))
    colors_for_df = [color_dict[x] for x in list(df[cell_type_key])]

    fig.axes[0].scatter(df[0], df[1],
                    s=50,
                    c=colors_for_df,
                    linewidth=1,
                    edgecolor='black', alpha=1)
    fig.savefig(f'{results_dir_path}/plots/umap.jpg',bbox_inches='tight')

# results quantification
def purity(tmp_ad, sample_name):
    SEACell_purity = SEACells.evaluate.compute_celltype_purity(tmp_ad, cell_type_key)
    plt.figure(figsize=(4,4))
    sns.boxplot(data=SEACell_purity, y='%s_purity' % cell_type_key)
    plt.title('Celltype Purity of SEACells computed from {} cells in {} ({})'.format(lineage_key, sample_name, study))
    sns.despine()
    plt.savefig(f'{results_dir_path}/plots/purity.jpg', bbox_inches='tight')
    return SEACell_purity

def compactness(tmp_ad, sample_name):
    compactness = SEACells.evaluate.compactness(tmp_ad, 'X_pca')
    plt.figure(figsize=(4,4))
    sns.boxplot(data=compactness, y='compactness')
    plt.title('Compactness of SEACells computed from {} cells in {} ({})'.format(lineage_key, sample_name, study))
    sns.despine()
    plt.savefig(f'{results_dir_path}/plots/compactness.jpg', bbox_inches='tight')
    return compactness

def separation(tmp_ad, sample_name):
    separation = SEACells.evaluate.separation(tmp_ad, 'X_pca',nth_nbr=1)
    plt.figure(figsize=(4,4))
    sns.boxplot(data=separation, y='separation')
    plt.title('Separation of SEACells computed from {} cells in {} ({})'.format(lineage_key, sample_name, study))
    sns.despine()
    plt.savefig(f'{results_dir_path}/plots/separation.jpg', bbox_inches='tight')
    return separation

def quantify_results(tmp_ad, sample_name):
    p = purity(tmp_ad, sample_name)
    c = compactness(tmp_ad, sample_name)
    s = separation(tmp_ad, sample_name)
    quant_df = pd.concat([p, c, s], axis=1)
    return quant_df

if __name__ == '__main__':
    # ----------------------------------------------------------------------
    # LOAD, FILTER, PREP
    # command line arguments
    ad_path = sys.argv[1]
    study = sys.argv[2]
    sample_name = sys.argv[3]
    results = sys.argv[4]
    split_path = sys.argv[5]

    results_dir_path = "{}/{}/{}".format(results, study, sample_name)
    try: os.mkdir("{}/{}".format(results, study))
    except Exception: pass
    try: os.mkdir(results_dir_path)
    except Exception: pass
    try: os.mkdir("%s/plots" % results_dir_path)
    except Exception: pass

    if len(sys.argv) > 6: lineage_key = sys.argv[6]
    else: lineage_key = "all"
    geneset_paths = ['/data/peer/sam/cytokine_central/references/gene_sets/spectra_dicts/spectra_genesetDict_LiningSublining_cytokineSignaling_v1_20230627.json']
    # geneset_paths = ['/Users/suhanibalachandran/Documents/mskwork/cytokine_central/references/gene_sets/spectra_dicts/spectra_genesetDict_MacMonoDC_cytokineSignaling_v3_20230623.json']
    # Column labels (change if they vary in this particular dataset)
    sample_key = 'sample'
    cell_type_key = 'subtype_broad'
    # Load anndata
    full_ad = sc.read(ad_path)
    split_ad = sc.read_mtx(split_path)
    full_ad.layers['original_X'] = full_ad.X
    full_ad.X = split_ad.X

# Subset to cell types of interest if necessary
    hematopoietic_cells = ['classical monocyte', 'macrophage', 'CD8-positive, alpha-beta T cell', 'alveolar macrophage', 'CD4-positive, alpha-beta T cell', 'CD1c-positive myeloid dendritic cell', 'regulatory T cell', 'mast cell', 'natural killer cell', 'B cell', 'plasma cell', 'dendritic cell','myeloid cell','conventional dendritic cell','non-classical monocyte','plasmacytoid dendritic cell','neutrophil']
    if lineage_key == 'hematopoietic':
        ad_subset = full_ad[ (full_ad.obs[sample_key] == sample_name) & (full_ad.obs[cell_type_key].isin(hematopoietic_cells))].copy()
    else: ad_subset = full_ad[ full_ad.obs[sample_key] == sample_name].copy()
        
    if len(ad_subset.obs_names) >= 500:
        # ----------------------------------------------------------------------
        # SEACELLS RUN
        try:
            # Filter genes
            sc.pp.filter_genes(ad_subset, min_counts=150)

            # Copy the counts to the "raw" attribute
            # raw_ad = sc.AnnData(ad_subset.X)
            # raw_ad.obs_names, raw_ad.var_names = ad_subset.obs_names, ad_subset.var_names
            # ad_subset.raw = raw_ad
            # del raw_ad
            ad_subset.X = ad_subset.layers['counts'].copy()

            # Normalize cells, log transform, compute highly variable genes
            sc.pp.normalize_total(ad_subset)
            sc.pp.log1p(ad_subset)
            sc.pp.highly_variable_genes(ad_subset, batch_key=sample_key, n_top_genes=2500, inplace=True)

            # Subset to highly variable genes + those of interest
            genes_list = ['ISG15', 'STAT1', 'NFKBIA', 'DRAM1', 'TREM2', 'SPP1', 'C1QA', 'CD163', 'SERPINE1', 'KLF10', 'HMOX1', 'MX1', 'CXCL9', 'IRF4', 'PTGS2', 'IL1B']
            conv_path = '/home/balachs2/work/references/symbol_to_ensembl_conversion_cytosig_20230724.csv'
            # conv_path = '/Users/suhanibalachandran/Documents/mskwork/cytokine_central/references/symbol_to_ensembl_conversion_cytosig_20230724.csv' # HARDCODED PATH - ensembl id conversion dictionary
            conv_df = pd.read_csv(conv_path)
            ref_anno = {}
            for p in geneset_paths:
                with open(p, 'rb') as infile: ref_anno.update(json.load(infile))
            genes_in_geneset = []
            for x in ref_anno.keys():
                for y in ref_anno[x].keys():
                    genes_in_geneset.extend(ref_anno[x][y])
            ensembl_ids = list(conv_df[conv_df['gene_symbol'].isin(genes_list + genes_in_geneset)]['ensembl_id'])
            genes_of_interest = ad_subset.var['highly_variable'].copy()
            for x in ensembl_ids:
                if x in genes_of_interest.index: genes_of_interest[x] = True

            # Compute principal components, neighbors, UMAP
            # THE BELOW ONLY WORKS FOR SCANPY V1.10
            # ad_subset.var['genes_of_interest'] = genes_of_interest
            # sc.pp.pca(ad_subset, mask_var='genes_of_interest')
            ad_subset.var['highly_variable_old'] = ad_subset.var['highly_variable'].copy()
            ad_subset.var['highly_variable'] = genes_of_interest
            sc.pp.pca(ad_subset, use_highly_variable=True)
            sc.pp.neighbors(ad_subset, use_rep='X_pca')
            sc.tl.umap(ad_subset)

            # Define model variables
            n_SEACells = max(ad_subset.shape[0]//75, 1)
            build_kernel_on = 'X_pca' # for scRNA-seq
            n_waypoint_eigs = min(n_SEACells, 10) # Number of eigenvalues to consider when initializing metacells
            waypoint_proportion = 0.9 # Proportion of metacells to initialize using waypoint analysis; currently not an available parameter

            # Run the model
            model = SEACells.core.SEACells(ad_subset,
                        build_kernel_on=build_kernel_on,
                        n_SEACells=n_SEACells,
                        n_waypoint_eigs=n_waypoint_eigs,
                        convergence_epsilon = 1e-5)

            # Initialize archetypes
            model.construct_kernel_matrix()
            model.initialize_archetypes()

            # Fit model
            model.fit(min_iter=10, max_iter=50)

            # Save summaries 
            a = model.get_hard_assignments()
            a.to_csv("{}/assignments.csv".format(results_dir_path))
            ad_subset.write('{}/{}_seacells.h5ad'.format(results_dir_path, sample_name))
            model.plot_convergence(save_as="{}/plots/convergence.jpg".format(results_dir_path), show=False)
            plot_umap(ad_subset, cell_type_key) # UMAP plot with SEACells over it

            # Quantify results (celltype purity, compactness, separation)
            quant_df = quantify_results(ad_subset, sample_name)

            # Summarize by (within-sample) SEACell
            summarized_ad = summarize_by_SEACell(ad_subset, SEACells_label='SEACell', summarize_layer='counts') # layer = X bc default is layer = raw which does not subset by var names
            summarized_ad.layers['raw'] = summarized_ad.X

            # Annotate summarized anndata with .obs information from the original single cell data
            summarized_ad.obs[sample_key] = [ad_subset.obs.loc[ad_subset.obs.index[ad_subset.obs['SEACell'] == x][0], sample_key] for x in summarized_ad.obs_names]
            summarized_ad.obs['# Single Cells'] = ad_subset.obs.groupby('SEACell').count().iloc[:,0].loc[summarized_ad.obs_names]
            summarized_ad.obs[cell_type_key] = ad_subset.obs.groupby('SEACell').apply(lambda x: x[cell_type_key].value_counts().index[0])
            summ_df = summarized_ad.obs.join(quant_df, rsuffix='_r')
            if (len(summ_df[cell_type_key].compare(summ_df['%s_r' % cell_type_key])) == 0):
                summ_df.drop('%s_r' % cell_type_key, axis=1, inplace=True)
            summ_df.reset_index(names='SEACell', inplace=True)
            summ_df.set_index([sample_key, 'SEACell'], inplace=True)
            summ_df.to_csv('{}/summary_statistics.csv'.format(results_dir_path))
            print("DONE!")
        except Exception as e:
            print(f'Sample {sample_name} failed with error {e}')
    else: print("DONE! less than 500 cells")



