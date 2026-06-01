import numpy as np
import pandas as pd
import scanpy as sc
import scanpy.external as sce
import os
os.environ["CUDA_VISIBLE_DEVICES"]="0"
import rapids_singlecell as rsc
import cupy as cp
# Enable `managed_memory`
import rmm
from rmm.allocators.cupy import rmm_cupy_allocator
rmm.reinitialize(
    managed_memory=True,
    pool_allocator=False,
)
cp.cuda.set_allocator(rmm_cupy_allocator)

sc.settings.verbosity = 0 # verbosity: errors (0), warnings (1), info (2), hints (3)
sc.settings.set_figure_params(dpi=100, 
                              fontsize=8, 
                              dpi_save=300, 
                              figsize=(5,5), 
                              format='pdf', 
                              vector_friendly=True, 
                              transparent=True)
sc.settings.figdir = '/path/to/your/fig/t'

sc.settings.n_jobs = 128
sc.settings.max_memory = 256
seed = 1234
np.random.seed(seed)

wd = '/path/to/your/wd'

adata_all = sc.read_h5ad(wd + '/all_cells.h5ad')

adata = adata_all[(adata_all.obs['anno_total'] == 'T cells')].copy()

adata.X = adata.layers['counts'].copy()

rsc.pp.filter_cells(adata, min_count = 200, qc_var='n_genes')
rsc.pp.filter_genes(adata, min_count = 20, qc_var='n_cells')

adata = adata[adata.obs.nFeature_RNA < 5000, :]
adata = adata[adata.obs.nCount_RNA < 15000, :]
adata = adata[adata.obs.percent_mt < 10, :]

# Set up un-normalized matrix to run scrublet on GPU
adata = adata.copy()
adata.X = adata.X.toarray()
rsc.get.anndata_to_GPU(adata)
rsc.pp.scrublet(adata, batch_key='sample', verbose = True)
print('\nscrublet results:', adata.obs['predicted_doublet'].value_counts())

# filter out doublets
adata = adata[adata.obs['predicted_doublet'] == False].copy()

# run normalization and log1p transformation on CPU
rsc.get.anndata_to_CPU(adata)
# Normalizing to median total counts
sc.pp.normalize_total(adata)
# Logarithmize the data:
print('\nrunning rsc.pp.log1p(adata)')
sc.pp.log1p(adata)

adata.layers['lognorm'] = adata.X.copy()

# Get the highly variable genes which are usually the most infromative genes.
sc.pp.highly_variable_genes(adata, 
                            n_top_genes=3000, 
                            batch_key="sample")

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
    adata.var.index.isin(total_gene),
    False, 
    adata.var['highly_variable']
)

# read cell cycle S genes and G2M genes
s_genes = pd.read_csv(wd + 'cell_cycle_s_genes.csv', header=0)['value'].tolist()

g2m_genes = pd.read_csv(wd + 'cell_cycle_g2m_genes.csv', header=0)['value'].tolist()

rsc.get.anndata_to_GPU(adata)

# Regress out effects of total counts per cell and the percentage of mitochondrial genes expressed.
rsc.pp.regress_out(adata, 
                   keys = ['total_counts_mt', 'nFeature_RNA', 
                           'nCount_RNA', 'percent_mt'])
# Scale the data to unit variance.
rsc.pp.scale(adata, max_value=10)

# PCA
rsc.pp.pca(adata, n_comps=50)
rsc.pp.harmony_integrate(adata, 
                         key=['sample', 'donor'], 
                         adjusted_basis = 'X_pca_harmony')

rsc.pp.neighbors(adata, 
                 n_neighbors=10, 
                 n_pcs=40, 
                 use_rep='X_pca_harmony', 
                 random_state=seed)

rsc.tl.leiden(adata, 
              resolution=0.6, 
              key_added='leiden_0_6', 
              random_state=seed)

rsc.get.anndata_to_CPU(adata)
sc.tl.paga(adata, 
           groups='leiden_0_6')
sc.pl.paga(adata, 
           plot=False, 
           random_state=seed)
sc.tl.umap(adata, 
           init_pos='paga', 
           random_state=seed)

adata.X = adata.layers['lognorm'].copy()

adata.write(wd + '/total_t_cells.h5ad')