from matplotlib import pyplot as plt
from pandas.errors import EmptyDataError
from scipy.stats import pearsonr
from numpy import linalg as LA
import networkx as nx
import seaborn as sns
import pandas as pd
import numpy as np
import random
import copy
import os

from util.centerality import plot_k_curve, right_target_global_centrality_t

plt.rcParams.update({'font.size': 17})


def getRanking(source_path, target_path, hormone = 'insulin', source_tissue='pancreas', target_tissue='skeletal_muscle'):
    """
    Read correlation matrix
    The correlation matrix represents spearman correlation between every pair of genes in source and target tissue.
    """
    print("\n[Reading correlation file ...]")
    corr_data = abs(pd.read_csv(f'data/{hormone}_eQTL_spearman_corr_matrix_with_snap_removed_CF_lt_R_latest_15k.csv'))
    #prob_data = pd.read_csv("../data/insulin_eQTL_spearman_prob_matrix_adj_removed_CF_lt_R_latest_10k.csv")
    
    print("\n[Correlation file loaded successfully.]")

    n = int(int(corr_data.shape[0])/2)

    genes = list(corr_data.columns[:n])
    
    print("\n[Raeding input files ...]")
    try:
        print(f'source path:{source_path}, {os.path.exists(source_path)}')
        print(f'target path:{target_path}, {os.path.exists(target_path)}')
        source_genes = list(pd.read_csv(source_path, header=None).iloc[:,0])
        target_genes = list(pd.read_csv(target_path, header=None).iloc[:,0])
    except EmptyDataError:
        print("\nFiles are empty or unable to read input csv files.")

    common_target_genes = np.intersect1d(genes, target_genes)
    common_source_genes = np.intersect1d(genes, source_genes)


    print("Number of target genes present in the data: ", len(common_target_genes))
    print("Number of source genes present in the data: ", len(common_source_genes))

    target_genes_indices = [i for i, e in enumerate(genes) if e in common_target_genes]
    source_genes_indices = [i for i, e in enumerate(genes) if e in common_source_genes]

    print("Number of target genes present in the data: ", len(target_genes_indices))
    print("Number of source genes present in the data: ", len(source_genes_indices))

    print(target_genes_indices)
    
    SNAP_data = pd.read_csv("data/SNAP_data/PPT-Ohmnet_gene_symbols.csv", index_col=0)
    #SNAP_data
    pancreas_df = SNAP_data.loc[SNAP_data['tissue'] == source_tissue].values
    skeletal_muscle_df = SNAP_data.loc[SNAP_data['tissue'] == target_tissue].values

    A_SNAP = np.zeros_like(corr_data.values, dtype=np.float32)

    #skeletal_muscle_adjacency = np.zeros_like((np.shape(skeletal_muscle_df)[0], 2))
    for i in range(np.shape(skeletal_muscle_df)[0]):
        n_offset = n
        try:
            from_gene = int(genes.index(skeletal_muscle_df[i,0])) + n_offset
            to_gene = int(genes.index(skeletal_muscle_df[i,1])) + n_offset    
            A_SNAP[from_gene, to_gene] = 1
            A_SNAP[to_gene, from_gene] = 1
        except:
            pass
            # print("In skeletal muscle ", skeletal_muscle_df[i, 0], " or ", skeletal_muscle_df[i,1], " not present.")             
        

    for i in range(np.shape(pancreas_df)[0]):
        n_offset = 0
        try:
            from_gene = int(genes.index(pancreas_df[i,0])) + n_offset
            to_gene = int(genes.index(pancreas_df[i,1])) + n_offset    
            A_SNAP[from_gene, to_gene] = 1
            A_SNAP[to_gene, from_gene] = 1
        except:
            pass
            # print("In pancreas ", pancreas_df[i, 0], " or ", pancreas_df[i,1], " not present.")  

    A = np.zeros_like(corr_data.values)
    A_orig = corr_data.values + A_SNAP
    A[:n,:n] = A_orig[n:,n:]
    A[:n,n:] = A_orig[n:,:n]
    A[n:,:n] = A_orig[:n,n:]
    A[n:,n:] = A_orig[:n,:n]

    print("Swapping over")
    l,g = right_target_global_centrality_t(A, num_layers=2, target_tissue = 1, target_gene_indices = source_genes_indices, p=0.9)
    # plt.hist(g[:n], bins = 'auto')

    # plot_df,results, filtered_results, lncRNA_results = plot_k_curve(genes, g, ground_truth_genes=common_target_genes, filtered=False, n=n)
    result = pd.DataFrame()
    result['Gene Name'] = genes
    result['Centrality'] = g[:n]
    
    return result
