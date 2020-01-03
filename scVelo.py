###### load python libraries
import sys
import scvelo as scv
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib as mpl
from matplotlib import cm
from matplotlib.colors import ListedColormap, LinearSegmentedColormap
from matplotlib import pyplot as plt
from matplotlib import rcParams

##### Setup scVelo
sc.settings.verbosity = 3
sc.set_figure_params(dpi=80, color_map='viridis')
scv.settings.set_figure_params('scvelo')


##### Load Data
LoomFile=sys.argv[1]
OUTFILE_PREFIX=sys.argv[2]

Experiment = OUTFILE_PREFIX
adata = scv.read(LoomFile, cache=True)
adata.var_names_make_unique()

cellnames = adata.obs_names
df = pd.DataFrame(adata.obs_names)
df.to_csv('Samples.txt', index=False)

anno = pd.read_csv('Annotations.txt')                ###################### Example of Annotations.txt is available on GitHub
adata.obs = anno


#### Basic Filtering
sc.pp.filter_cells(adata, min_genes=0)
sc.pp.filter_genes(adata, min_cells=0)

##### Normalization
sc.pp.normalize_per_cell(adata, counts_per_cell_after=1e4)
sc.pp.log1p(adata)
adata.raw = adata


####### PCA, Clustering and UMAP

########## set n_neighbors, n_pcs and n_comps based on the number of samples#######


sc.tl.pca(adata, svd_solver='arpack', n_comps=3)
sc.pp.neighbors(adata, n_neighbors=2, n_pcs=3)                 
sc.tl.umap(adata)
sc.tl.louvain(adata)
FigureOut = Experiment + "_UMAP_scVelo.png"
rcParams['figure.figsize'] = 15, 15
sc.pl.umap(adata, color=['louvain'], save=FigureOut)
FigureOut = Experiment + "_UMAP_SampleType_scVelo.png"
rcParams['figure.figsize'] = 15, 15
sc.pl.umap(adata, color=['SampleType'], save=FigureOut)

###### scVelo Velocity Estimation
scv.pp.moments(adata, n_pcs=3)
scv.tl.velocity(adata, mode='stochastic')
scv.tl.velocity_graph(adata)
scv.tl.velocity_embedding(adata, basis='umap')

##### scVelo Visualization
FigureOut = Experiment + "_VelocityEmbedding_scVelo.png"
rcParams['figure.figsize'] = 15, 15
scv.pl.velocity_embedding(adata, basis='umap', save=FigureOut)

FigureOut = Experiment + "_VelocityStream_scVelo.png"
rcParams['figure.figsize'] = 15, 15
scv.pl.velocity_embedding_stream(adata, basis='umap', save=FigureOut)

#### Variance Velocity of Genes of Interest
FigureOut = Experiment + "_Cited1-Wnt4-Hnf4a_Velocity_scVelo.png"
scv.pl.velocity(adata, var_names=['Notch2', 'Wnt4', 'Hnf4a'])            ###### Change genes based on your interest

#### writing out scVelo object
results_file = Experiment + 'scVeloAnalysis.h5ad'
adata.write(results_file, compression='gzip')