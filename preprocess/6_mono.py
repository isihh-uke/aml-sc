import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce
import os
#os.environ["CUDA_VISIBLE_DEVICES"]="1"

#import torch
#torch.cuda.is_available()

#import rapids_singlecell as rsc
from matplotlib.pyplot import rc_context
import matplotlib.pyplot as plt
import seaborn as sns
import math
import sys

# Enable `pool_allocator`
#import cupy as cp
# Enable `managed_memory`
#import rmm
#from rmm.allocators.cupy import rmm_cupy_allocator
#rmm.reinitialize(
#    managed_memory=True,
#    pool_allocator=False,
#)
#cp.cuda.set_allocator(rmm_cupy_allocator)

# Enable `pool_allocator`
#import rmm
#from rmm.allocators.cupy import rmm_cupy_allocator
#rmm.reinitialize(
#    managed_memory=False,
#    pool_allocator=True,
#)
#cp.cuda.set_allocator(rmm_cupy_allocator)

sc.settings.verbosity = 0 # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.set_figure_params(dpi=100, fontsize=8, dpi_save=300, figsize=(5,5), format='pdf', vector_friendly = True, transparent=True)
sc.settings.figdir = '/path/to/your/fig/mono/'

sc.settings.n_jobs = 128
sc.settings.max_memory = 256

seed = 1234
np.random.seed(seed)

wd = '/path/to/your/wd'

mono = ['cMono', 'ncMono']

adata_my = sc.read_h5ad(wd + 'myeloid_cells.h5ad')

adata_mo = adata_my[adata_my.obs['anno_my'].isin(mono)].copy()

adata_mo

adata = adata_mo.copy()
#rsc.get.anndata_to_GPU(adata)
# get the highly variable genes which are usually the most infromative genes.
sc.pp.highly_variable_genes(adata, n_top_genes=3000, batch_key="sample")

# Get the list of all variable names (gene names from adata.var_names) 
all_gene_names = adata.var_names.tolist()

# hemoglobin genes
hb_gene = adata.var[adata.var.highly_variable == True].index[adata.var[adata.var.highly_variable == True].index.str.contains(("^HB[^(P)]"))].tolist()
# un-annotated genes
ens_gene = adata.var[adata.var.highly_variable == True].index[adata.var[adata.var.highly_variable == True].index.str.startswith(('ENS'))].tolist()
# mito genes
mt_gene = adata.var[adata.var.highly_variable == True].index[adata.var[adata.var.highly_variable == True].index.str.startswith(('MT-'))].tolist()
# hsp genes
hsp_gene = adata.var[adata.var.highly_variable == True].index[adata.var[adata.var.highly_variable == True].index.str.startswith(('HSP'))].tolist()
# ribosomal genes
rp_gene = adata.var[adata.var.highly_variable == True].index[adata.var[adata.var.highly_variable == True].index.str.startswith(('RPS', 'RPL'))].tolist()

total_gene = hb_gene + ens_gene + mt_gene + hsp_gene + rp_gene

# Use np.where to update the 'highly_variable' column in adata.var
adata.var['highly_variable'] = np.where(
    adata.var.index.isin(total_gene),  # Condition: Check if the gene names (index of adata.var) are in total_gene
    False,                              # Value if condition is True: Set to False
    adata.var['highly_variable']       # Value if condition is False: Keep the current value of 'highly_variable'
)

# read cell cycle S genes and G2M genes
s_genes = pd.read_csv(wd+'cell_cycle_s_genes.csv', header=0)['value'].tolist()

g2m_genes = pd.read_csv(wd+'cell_cycle_g2m_genes.csv', header=0)['value'].tolist()

sc.tl.score_genes(adata, gene_list = s_genes, score_name = 'S_score', random_state = seed)
sc.tl.score_genes(adata, gene_list = g2m_genes, score_name = 'G2M_score', random_state = seed)

print('adata size before scaling: ', sys.getsizeof(adata)/(pow(1024, 3)), 'GB') 

# Regress out effects of total counts per cell and the percentage of mitochondrial genes expressed. Scale the data to unit variance.
sc.pp.regress_out(adata, keys = ['total_counts', 'total_counts_mt', 'nFeature_RNA', 'nCount_RNA', 'percent_mt', 'S_score', 'G2M_score'])
sc.pp.scale(adata, max_value=10)
print('adata size after scaling: ', sys.getsizeof(adata)/(pow(1024, 3)), 'GB') 

# pca
sc.pp.pca(adata, n_comps=50)
sc.external.pp.harmony_integrate(adata, key=['sample', 'donor'], adjusted_basis = 'X_pca_harmony')
#adata = adata.raw.to_adata()
print('adata size after harmony: ', sys.getsizeof(adata)/(pow(1024, 3)), 'GB') 

#seed = 1234
#np.random.seed(seed)
sc.pp.neighbors(adata, n_neighbors=10, n_pcs=40, use_rep='X_pca_harmony', random_state=seed)
sc.tl.leiden(adata, resolution=0.6, key_added='leiden_0_6', random_state=seed)
sc.tl.paga(adata, groups='leiden_0_6')
sc.pl.paga(adata, plot=False, random_state=seed)  # remove `plot=False` if you want to see the coarse-grained graph
sc.tl.umap(adata, init_pos='paga', random_state=seed)
print('adata size after leiden: ', sys.getsizeof(adata)/(pow(1024, 3)), 'GB') 

adata.X = adata.layers['lognorm'].copy()

print('adata size after removing scaled mtx: ', sys.getsizeof(adata)/(pow(1024, 3)), 'GB') 


adata.write(wd + 'all_cells_h5ad/mono.h5ad')