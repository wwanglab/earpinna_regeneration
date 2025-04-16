# -*- coding: utf-8 -*-
'''
Description: 
'''

import os,sys
import scanpy as sc
import numpy as np
import matplotlib.pyplot as plt
import anndata
import random
import pandas as pd
from scvi.external import GIMVI
import scvi

def getDefaultColors(n, type = 1):
    if type == 1:
        colors = ["#ff1a1a", "#1aff1a", "#1a1aff", "#ffff1a", "#ff1aff",
                "#ff8d1a", "#7cd5c8", "#c49a3f", "#5d8d9c", "#90353b",
                "#507d41", "#502e71", "#1B9E77", "#c5383c", "#0081d1",
                "#674c2a", "#c8b693", "#aed688", "#f6a97a", "#c6a5cc",
                "#798234", "#6b42c8", "#cf4c8b", "#666666", "#ffd900",
                "#feb308", "#cb7c77", "#68d359", "#6a7dc9", "#c9d73d"]
    return colors


args = sys.argv
spdata = args[1]
scdata = args[2]
mask = args[3]
outdir = args[4]
os.system("mkdir -p %s"%(outdir))
os.chdir(outdir)
st_adata = sc.read_h5ad(spdata)
scvi.settings.num_threads = 30
sc_adata = sc.read_h5ad(scdata)
sc_adata.__dict__['_raw'].__dict__['_var'] = sc_adata.__dict__['_raw'].__dict__['_var'].rename(columns={'_index': 'features'})

sc.pp.filter_cells(sc_adata, min_counts = 10)
#only use genes in both datasets
G = 3000

sc_adata.layers["counts"] = sc_adata.raw.X.astype(int).copy()
sc.pp.filter_cells(st_adata, min_counts = 10)

st_adata.layers["counts"] = st_adata.X.astype(int).copy()
sc.pp.filter_cells(st_adata, min_counts = 10)


# filter genes to be the same on the spatial data
intersect = np.intersect1d(sc_adata.var_names, st_adata.var_names)
st_adata = st_adata[:, intersect].copy()
sc_adata = sc_adata[:, intersect].copy()
G = len(intersect)

GIMVI.setup_anndata(st_adata,layer="counts")
GIMVI.setup_anndata(sc_adata, layer="counts", labels_key="celltype",batch_key="time")
model = GIMVI(sc_adata, st_adata)
model.train(200)
GIMVI.save(model, dir_path=outdir, overwrite=True, save_anndata=True)
model = GIMVI.load(outdir)
st_adata = sc.read_h5ad(outdir+'/adata_spatial.h5ad')
sc_adata = sc.read_h5ad(outdir+'/adata_seq.h5ad')
_, imputed = model.get_imputed_values(normalized=True)
st_adata.X = imputed
sc.pp.log1p(st_adata)
sc.pp.highly_variable_genes(st_adata, flavor="seurat", n_top_genes=2000)
sc.pp.pca(st_adata)
sc.pp.neighbors(st_adata)
sc.tl.umap(st_adata)
sc.tl.leiden(st_adata, key_added="res.0.4",resolution=0.4)
sc.tl.leiden(st_adata, key_added="res.0.8",resolution=0.8)
sc.tl.leiden(st_adata, key_added="res.1",resolution=1)
sc.tl.leiden(st_adata, key_added="res.1.4",resolution=1.4)
st_adata.write(outdir+"/cellbin_imputated.h5ad")
plt.rcParams["figure.figsize"] = (8, 8)
fig, axs = plt.subplots(2, 2, figsize=(16, 16))
sc.pl.umap(st_adata, color="res.0.4", wspace=0.4,show=False,ax=axs[0][0])
sc.pl.umap(st_adata, color="res.0.8", wspace=0.4,show=False,ax=axs[0][1])
sc.pl.umap(st_adata, color="res.1", wspace=0.4,show=False,ax=axs[1][0])
sc.pl.umap(st_adata, color="res.1.4", wspace=0.4,show=False,ax=axs[1][1])
plt.savefig(outdir+"/cluster_umap.pdf",bbox_inches="tight")
fig, axs = plt.subplots(2, 2, figsize=(25, 16))
cluster_number = st_adata.obs['res.0.4'].cat.categories.shape[0]
colors = getDefaultColors(cluster_number, type = 2)
sc.pl.spatial(st_adata, img_key="hires", color="res.0.4", size=1.5,spot_size = 20,palette=colors,show=False,ax=axs[0][0])
cluster_number = st_adata.obs['res.0.8'].cat.categories.shape[0]
colors = getDefaultColors(cluster_number, type = 2)
sc.pl.spatial(st_adata, img_key="hires", color="res.0.8", size=1.5,spot_size = 20,palette=colors,show=False,ax=axs[0][1])
cluster_number = st_adata.obs['res.1'].cat.categories.shape[0]
colors = getDefaultColors(cluster_number, type = 2)
sc.pl.spatial(st_adata, img_key="hires", color="res.1", size=1.5,spot_size = 20,palette=colors,show=False,ax=axs[1][0])
cluster_number = st_adata.obs['res.1.4'].cat.categories.shape[0]
colors = getDefaultColors(cluster_number, type = 2)
sc.pl.spatial(st_adata, img_key="hires", color="res.1.4", size=1.5,spot_size = 20,palette=colors,show=False,ax=axs[1][1])
plt.savefig("spatial_cluster_umap.pdf",bbox_inches="tight")
resolution = 'res.1' #choose best resolution
res = pd.DataFrame(st_adata.obs, columns = ["x", "y",resolution], index = st_adata.obs.index)
res.to_csv("bin1clu.txt",sep = '\t',index =False)
clusters = st_adata.obs[resolution].cat.categories
cluster_number = clusters.shape[0]
colors = getDefaultColors(cluster_number, type = 2)
flout = open("color.list",'w')
for i in range(cluster_number):
    if clusters[i] == 'low_quality':
        flout.write(clusters[i] + '\t#ffffff\n')
    else:
        flout.write(clusters[i] + '\t' + colors[i] + '\n')
flout.close()
os.system('cell_bin_plot bin1clu.txt %s color.list cluster_plot.tif'%(mask))
os.system('mkdir -p cluster_split')
os.chdir("cluster_split")
flout = open("color.list",'w')
flout.write('1\t#ff0000\n')
flout.write('low_quality\t#ffffff\n')
flout.close()
for i in clusters:
    if i == 'low_quality':
        continue
    tmp = res.copy()
    tmp.loc[tmp[resolution] != i, [resolution]] = 'low_quality'
    tmp.loc[tmp[resolution] == i, [resolution]] = '1'
    tmp.to_csv("bin1clu_%s.txt"%(i),sep = '\t',index =False)
    os.system('cell_bin_plot bin1clu_%s.txt %s color.list cluster_plot_%s.tif'%(i, mask, i))
#st_adata = removeBiasGenes(st_adata)
if st_adata.raw:
    del st_adata.raw
if 'log1p' in st_adata.uns.keys():
    st_adata.uns['log1p']["base"] = None
sc.pp.highly_variable_genes(st_adata, flavor="seurat", n_top_genes=2000)
sc.tl.rank_genes_groups(st_adata, resolution, method='wilcoxon')
plt.rcParams["figure.figsize"] = (7,7)
sc.pl.rank_genes_groups(st_adata, n_genes=25, sharey=False)
plt.savefig("rank_genes_groups.pdf")

markers = pd.DataFrame(st_adata.uns['rank_genes_groups']['names']).head(3).stack().values.tolist()
markers = list(set(markers))
sc.pl.stacked_violin(st_adata, markers, groupby=resolution, rotation=90)
plt.savefig("marker_genes_violin.pdf")

degs = pd.DataFrame(st_adata.uns['rank_genes_groups']['names']).head(50)
degs.to_csv("degs.txt",sep = '\t')



