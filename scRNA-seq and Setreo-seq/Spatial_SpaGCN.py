###########################################################################
######################2. Import python modules#############################
###########################################################################

import os,csv,re,sys
import pandas as pd
import numpy as np
import scanpy as sc
import math
import SpaGCN as spg
from scipy.sparse import issparse
from scipy.spatial import distance#zuolulu
import random, torch
import warnings
warnings.filterwarnings("ignore")
import matplotlib.colors as clr
import matplotlib.pyplot as plt
import SpaGCN as spg
import cv2

args = sys.argv
indir=args[1]
pre=args[2]
outdir=args[3]
clunum=int(args[4])

###########################################################################
##############################3. Read in data##############################
###########################################################################
adata=sc.read(indir+"/"+pre+"_imputated.h5ad")

###########################################################################
#####4. Integrate gene expression and histology into a Graph###############
###########################################################################
#Set coordinates
x_array = adata.obsm['spatial'][:,0]
y_array = adata.obsm['spatial'][:,1]
X = np.array([x_array, y_array]).T.astype(np.float32)
adj = distance.cdist(X, X, 'euclidean') #adjmatrix
np.savetxt(outdir+"/"+pre+'_imputated_adj.csv', adj, delimiter=',')

###########################################################################
###############5.Spatial domain detection using SpaGCN#####################
###########################################################################
###############################5.1 Expression data preprocessing
adata=sc.read(indir+"/"+pre+"_imputated.h5ad")
adj=np.loadtxt(outdir+"/"+pre+'_imputated_adj.csv', delimiter=',')
adata.var_names_make_unique()
spg.prefilter_genes(adata,min_cells=3) # avoiding all genes are zeros
spg.prefilter_specialgenes(adata)

###############################5.2 Set hyper-parameters
#Set hyper-parameters
p=0.5 #Percentage of total expression contributed by neighborhoods
#Find the l value given p
l=spg.search_l(p, adj, start=0.01, end=1000, tol=0.01, max_run=100)

n_clusters=clunum
r_seed=t_seed=n_seed=100
res=spg.search_res(adata, adj, l, n_clusters, start=0.7, step=0.1, tol=5e-3, lr=0.001, max_epochs=20, r_seed=r_seed, t_seed=t_seed, n_seed=n_seed)

################################5.3 Run SpaGCN
clf=spg.SpaGCN()
clf.set_l(l)
#Set seed
random.seed(r_seed)
torch.manual_seed(t_seed)
np.random.seed(n_seed)
#Run
clf.train(adata,adj,init_spa=True,init="louvain",res=res, tol=5e-3, lr=0.001, max_epochs=200)
y_pred, prob=clf.predict()
adata.obs["pred"]= y_pred
adata.obs["pred"]=adata.obs["pred"].astype('category')
#Do cluster refinement(optional)
#shape="hexagon" for Visium data, "square" for ST data.
adj_2d=spg.calculate_adj_matrix(x=x_array,y=y_array, histology=False)
refined_pred=spg.refine(sample_id=adata.obs.index.tolist(), pred=adata.obs["pred"].tolist(), dis=adj_2d, shape="hexagon")
adata.obs["refined_pred"]=refined_pred
adata.obs["refined_pred"]=adata.obs["refined_pred"].astype('category')
#Save results
adata.write_h5ad(outdir+"/"+pre+"_imputated_SpaGCN.h5ad")

################################5.4 Plot spatial domains
adata=sc.read(outdir+"/"+pre+"_imputated_SpaGCN.h5ad")
#Set colors used
plot_color=["#F56867","#FEB915","#C798EE","#59BE86","#7495D3","#D1D1D1","#6D1A9C","#15821E","#3A84E6","#997273","#787878","#DB4C6C","#9E7A7A","#554236","#AF5F3C","#93796C","#F9BD3F","#DAB370","#877F6C","#268785"]
#Plot spatial domains
domains="pred"
num_celltype=len(adata.obs[domains].unique())
adata.uns[domains+"_colors"]=list(plot_color[:num_celltype])
ax=sc.pl.scatter(adata,alpha=1,x="imagerow",y="imagecol",color=domains,title=domains,color_map=plot_color,show=False,size=100000/adata.shape[0])
ax.set_aspect('equal', 'box')
ax.axes.invert_yaxis()
plt.savefig(outdir+"/"+pre+"_imputated_SpaGCN_pred.png", dpi=600)
plt.close()

#Plot refined spatial domains
domains="refined_pred"
num_celltype=len(adata.obs[domains].unique())
adata.uns[domains+"_colors"]=list(plot_color[:num_celltype])
ax=sc.pl.scatter(adata,alpha=1,x="imagerow",y="imagecol",color=domains,title=domains,color_map=plot_color,show=False,size=100000/adata.shape[0])
ax.set_aspect('equal', 'box')
ax.axes.invert_yaxis()
plt.savefig(outdir+"/"+pre+"_imputated_SpaGCN_refined_pred.png", dpi=600)
plt.close()
