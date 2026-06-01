import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce
import os
os.environ["CUDA_VISIBLE_DEVICES"]="1"

#import torch
#torch.cuda.is_available()

import rapids_singlecell as rsc
from matplotlib.pyplot import rc_context
import matplotlib.pyplot as plt
import seaborn as sns
import math
import sys

# Enable `pool_allocator`
import cupy as cp
# Enable `managed_memory`
import rmm
from rmm.allocators.cupy import rmm_cupy_allocator
rmm.reinitialize(
    managed_memory=True,
    pool_allocator=False,
)
cp.cuda.set_allocator(rmm_cupy_allocator)

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
sc.settings.figdir = '/home/big/zheng_song/aml_python/com_s15_rerun/figs/harmony/hsc_prolif'

sc.settings.n_jobs = 128
sc.settings.max_memory = 256

seed = 1234
np.random.seed(seed)

wd = '/home/big/zheng_song/aml_python/com_s15_rerun/'

# read all cells adata
adata_all = sc.read_h5ad(wd + 'all_cells_h5ad/all_cells_id_corrected_meta_data_raw_counts_harmony.h5ad')
# subset hsc prolif part
adata_hsc = adata_all[adata_all.obs['anno_total'].isin(['Prolif. cells', 'HSCs'])]

adata = adata_hsc.copy()

del adata_all
del adata_hsc

# copy raw counst into X 
adata.X = adata.layers['counts'].copy()
print('adata before filtering:', adata)

with rc_context({'figure.figsize': (4, 3)}):
    ax = sc.pl.scatter(adata, x = 'nCount_RNA', y = 'nFeature_RNA', color = 'percent_mt', 
                       title='', size=5, show = False)
    _ = ax.set_xlabel(xlabel= 'Number of UMIs', fontsize = 10)
    _ = ax.set_ylabel(ylabel= 'Number of genes', fontsize = 10)
    _ = ax.set_title(label="% of mitochondrial genes")
    _ = ax.grid(False)
    plt.savefig(str(sc.settings.figdir) + '/counts_genes_mito_hsc_prolif_before_qc.pdf', dpi=300, bbox_inches='tight')

with rc_context({'figure.figsize': (8, 2.5), 'font.size':10}):
    ax = sc.pl.violin(adata, ['nFeature_RNA'], groupby='sample', 
                 log=False, jitter=True, scale='width',
                 inner = 'box', legend = False, dodge = False, fill=False, show=False)
    _ = ax.grid(False)
    _ = ax.set_xlabel('')
    _ = ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
    _ = ax.set_ylabel('Number of genes')
    _ = ax.set_yticklabels(ax.get_yticklabels())
    plt.savefig(str(sc.settings.figdir) + '/hsc_prolif_feature_violin_before_qc.pdf', dpi=300, bbox_inches='tight')

with rc_context({'figure.figsize': (8, 2.5), 'font.size':10}):
    ax = sc.pl.violin(adata, ['percent_mt'], groupby='sample', 
                 log=False, jitter=True, scale='width',
                 inner = 'box', legend = False, dodge = False, fill=False, show=False)
    _ = ax.grid(False)
    _ = ax.set_xlabel('')
    _ = ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
    _ = ax.set_ylabel('% of mitochondrial genes')
    _ = ax.set_yticklabels(ax.get_yticklabels())
    plt.savefig(str(sc.settings.figdir) + '/hsc_prolif_mt_violin_before_qc.pdf', dpi=300, bbox_inches='tight')

with rc_context({'figure.figsize': (8, 2.5), 'font.size':10}):
    ax = sc.pl.violin(adata, ['nCount_RNA'], groupby='sample', 
                 log=False, jitter=True, scale='width',
                 inner = 'box', legend = False, dodge = False, fill=False, show=False)
    _ = ax.grid(False)
    _ = ax.set_xlabel('')
    _ = ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
    _ = ax.set_ylabel('Number of UMIs')
    _ = ax.set_yticklabels(ax.get_yticklabels())
    plt.savefig(str(sc.settings.figdir) + '/hsc_prolif_count_violin_before_qc.pdf', dpi=300, bbox_inches='tight')

# filter low frequence genes and low counts cells
sc.pp.filter_cells(adata, min_genes = 200)
sc.pp.filter_genes(adata, min_cells = 20)

# filter our low quality cells
adata = adata[adata.obs.nFeature_RNA < 5000, :]
adata = adata[adata.obs.nCount_RNA < 15000, :]
adata = adata[adata.obs.percent_mt < 10, :]

with rc_context({'figure.figsize': (4, 3)}):
    ax = sc.pl.scatter(adata, x = 'nCount_RNA', y = 'nFeature_RNA', color = 'percent_mt', 
                       title='', size=5, show = False)
    _ = ax.set_xlabel(xlabel= 'Number of UMIs', fontsize = 10)
    _ = ax.set_ylabel(ylabel= 'Number of genes', fontsize = 10)
    _ = ax.set_title(label="% of mitochondrial genes")
    _ = ax.grid(False)
    plt.savefig(str(sc.settings.figdir) + '/counts_genes_mito_hsc_prolif_after_qc.pdf', dpi=300, bbox_inches='tight')

with rc_context({'figure.figsize': (8, 2.5), 'font.size':10}):
    ax = sc.pl.violin(adata, ['nFeature_RNA'], groupby='sample', 
                 log=False, jitter=True, scale='width',
                 inner = 'box', legend = False, dodge = False, fill=False, show=False)
    _ = ax.grid(False)
    _ = ax.set_xlabel('')
    _ = ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
    _ = ax.set_ylabel('Number of genes')
    _ = ax.set_yticklabels(ax.get_yticklabels())
    plt.savefig(str(sc.settings.figdir) + '/hsc_prolif_feature_violin_after_qc.pdf', dpi=300, bbox_inches='tight')

with rc_context({'figure.figsize': (8, 2.5), 'font.size':10}):
    ax = sc.pl.violin(adata, ['percent_mt'], groupby='sample', 
                 log=False, jitter=True, scale='width',
                 inner = 'box', legend = False, dodge = False, fill=False, show=False)
    _ = ax.grid(False)
    _ = ax.set_xlabel('')
    _ = ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
    _ = ax.set_ylabel('% of mitochondrial genes')
    _ = ax.set_yticklabels(ax.get_yticklabels())
    plt.savefig(str(sc.settings.figdir) + '/hsc_prolif_mt_violin_after_qc.pdf', dpi=300, bbox_inches='tight')

with rc_context({'figure.figsize': (8, 2.5), 'font.size':10}):
    ax = sc.pl.violin(adata, ['nCount_RNA'], groupby='sample', 
                 log=False, jitter=True, scale='width',
                 inner = 'box', legend = False, dodge = False, fill=False, show=False)
    _ = ax.grid(False)
    _ = ax.set_xlabel('')
    _ = ax.set_xticklabels(ax.get_xticklabels(), rotation=45, ha='right')
    _ = ax.set_ylabel('Number of UMIs')
    _ = ax.set_yticklabels(ax.get_yticklabels())
    plt.savefig(str(sc.settings.figdir) + '/hsc_prolif_count_violin_after_qc.pdf', dpi=300, bbox_inches='tight')

print('\nadata after filtering:\n\n', adata)

adata = adata.copy()

rsc.get.anndata_to_GPU(adata)
#adata.write(wd + 'all_cells_h5ad/hsc_prolif_cells_id_corrected_meta_data_raw_counts_harmony_re_clean.h5ad')

rsc.pp.scrublet(adata, batch_key='sample', verbose = True)

print('\nscrublet results:', adata.obs['predicted_doublet'].value_counts())

# we filter out doublets
adata = adata[adata.obs['predicted_doublet'] == False].copy()

print('\nadata after removing dbls:\n\n', adata)

#rsc.get.anndata_to_CPU(adata)

# Normalizing to median total counts
rsc.pp.normalize_total(adata)
# Logarithmize the data:
rsc.pp.log1p(adata)

adata.layers['lognorm'] = adata.X.copy()


# get the highly variable genes which are usually the most infromative genes.
rsc.pp.highly_variable_genes(adata, n_top_genes=3000, batch_key="sample")

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
rsc.pp.regress_out(adata, keys = ['total_counts', 'total_counts_mt', 'nFeature_RNA', 'nCount_RNA', 'percent_mt'])

# scale matrix for dimentional reduction
rsc.pp.scale(adata, max_value=10)
print('adata size after scaling: ', sys.getsizeof(adata)/(pow(1024, 3)), 'GB') 

#adata.write(wd + 'all_cells_h5ad/hsc_prolif_id_corrected_meta_data_raw_counts_harmony_re_clean_scaled.h5ad')

# pca
rsc.pp.pca(adata, n_comps=50)
rsc.pp.harmony_integrate(adata, key=['sample', 'donor'], adjusted_basis = 'X_pca_harmony')
#adata = adata.raw.to_adata()
adata.write(wd + 'all_cells_h5ad/hsc_prolif_id_corrected_meta_data_raw_counts_harmony_re_clean_re_harmony.h5ad')
print('adata size after harmony: ', sys.getsizeof(adata)/(pow(1024, 3)), 'GB') 

#seed = 1234
#np.random.seed(seed)
rsc.pp.neighbors(adata, n_neighbors=10, n_pcs=20, use_rep='X_pca_harmony', random_state=seed)
rsc.tl.leiden(adata, resolution=0.6, key_added='leiden_0_6', random_state=seed)
rsc.get.anndata_to_CPU(adata)
sc.tl.paga(adata, groups='leiden_0_6')
sc.pl.paga(adata, plot=False, random_state=seed)  # remove `plot=False` if you want to see the coarse-grained graph
sc.tl.umap(adata, init_pos='paga', random_state=seed)
print('adata size after leiden: ', sys.getsizeof(adata)/(pow(1024, 3)), 'GB') 

adata.X = adata.layers['lognorm'].copy()
                      
adata.write(wd + 'all_cells_h5ad/hsc_prolif_id_corrected_meta_data_raw_counts_harmony_re_clean_re_harmony_leiden.h5ad')
                      
