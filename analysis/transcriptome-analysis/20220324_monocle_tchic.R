# To start: load the necessary packages
#install.packages("Matrix")
library(dplyr)
library(Seurat)
library(patchwork)
library(monocle3)
library(ggplot2)
library(dplyr)
library(Matrix)

# Load your dataset
pathToData = "/Users/m.blotenburg/Documents/Projects/TCHIC/data/jupyter_notebooks/rep234-seqrun5651-6347/"
umap_coord<- read.csv(file= paste0(pathToData,"umap_coord.csv"),sep=",", header = TRUE, row.names = 1)
adata_obs <- read.csv(file= paste0(pathToData,"adata_obs.csv"),sep=",", header = TRUE, row.names = 1)
adata_var <- read.csv(file= paste0(pathToData,"adata_var.csv"),sep=",", header = TRUE, row.names = 1)
head(umap_coord)
head(adata_obs)
head(adata_var)
raw_counts <- readMM(paste0(pathToData,"sparse_matrix_int.mtx"))
head(raw_counts)

########

# transfer labels
row.names(umap_coord) <- row.names(adata_obs)
adata_var['gene_short_name'] <- row.names(adata_var)

############ MONOCLE ############
# import cds
cds <- new_cell_data_set(raw_counts,
                         cell_metadata = adata_obs,
                         gene_metadata = adata_var)

cds
## Step 1: Normalize and pre-process the data
cds <- preprocess_cds(cds, num_dim = 100)

# check PCA
plot_pc_variance_explained(cds)

## Step 2: Reduce the dimensions using UMAP
# you can use tsne by changing reduction_method to "tSNE"
cds <- reduce_dimension(cds, reduction_method = "UMAP",
                        preprocess_method = "PCA")

# check umap
# trajectory is not yet calculated: show_trajectory_graph = F
# cell_size defines dot size (default is very small: 0.35) - same for group_label_size
# color_cells_by can be any interesting feature
plot_cells(cds, show_trajectory_graph = F,
           cell_size = 1, color_cells_by = "celltype",
           group_label_size = 4)

plot_cells(cds, show_trajectory_graph = F,
           cell_size = 1, color_cells_by = "day",
           group_label_size = 4)

# or plot expression of some random genes:
plot_cells(cds, show_trajectory_graph = F,
           genes=c("Sox2"
                   ,"Sox17"
                   ,"Gata6"
                   , "Hand2"),
           cell_size = 1)

## optional: batch effect removal
# It is a good idea to check your plate distribution with 'plot_cells' to see if there is a plate bias (= batch effect)
# define your batch in 'alignment_group'
## Remove batch effects with cell alignment
# I skipped this
# cds <- align_cds(cds, preprocess_method = "PCA", alignment_group = "batch")


## Step 3: Cluster the cells
# cluster_method can be either "louvain" or "leiden"
# i set the k lower since the amount of cells is low. default is 20
# you can play around with the k (# nearest neighbours used for clustering) and see what happens with your clusters
cds <- cluster_cells(cds, reduction_method = "UMAP",
                     cluster_method = "leiden",
                     k = 10)

# check your new clusters
plot_cells(cds, show_trajectory_graph = F, color_cells_by = "cluster",
           cell_size = 1,
           group_label_size = 4)

plot_cells(cds, show_trajectory_graph = F,
           cell_size = 1, color_cells_by = "celltype",
           group_label_size = 4)

# compare to seurat
seurat <- plot_cells(cds, show_trajectory_graph = F,
                     cell_size = 1, color_cells_by = "celltype",
                     group_label_size = 4)
monocle <- plot_cells(cds, show_trajectory_graph = F,
                      cell_size = 1,
                      group_label_size = 4)
seurat + monocle

# transfer back original UMAP coordinates
row.names(umap_coord) <- row.names(adata_obs)
cds@int_colData@listData$reducedDims$UMAP <- umap_coord

#cds@int_colData@listData$reducedDims$UMAP

# check UMAP
plot_cells(cds, show_trajectory_graph = F,
           cell_size = 1, color_cells_by = "celltype",
           group_label_size = 4)

## Step 4: Learn a graph
cds <- learn_graph(cds, use_partition = T)

plot_cells(cds,
           color_cells_by = "celltype",cell_size=1,
           label_groups_by_cluster=F,
           label_leaves=F,
           label_branch_points=F)

## Step 5: Order cells in pseudotime - define root node manually
cds <- order_cells(cds)

# plot pseudotime
plot_cells(cds,
           color_cells_by = "pseudotime", cell_size = 1,
           label_cell_groups=FALSE,
           label_leaves=FALSE,
           label_branch_points=FALSE,
           graph_label_size=1.5)


colnames(cds@colData)
## plots
# check if clustering is biased by number of counts
plot_cells(cds, color_cells_by = "n_counts",
           show_trajectory_graph = F,
           cell_size = 1)

# check if clustering is biased by number of genes covered
plot_cells(cds, color_cells_by = "n_genes",
           show_trajectory_graph = F,
           cell_size = 1)

# check if clustering is biased by mitochondrial reads
plot_cells(cds, color_cells_by = "percent_mito",
           show_trajectory_graph = F,
           cell_size = 1)

# Check seurat clusters
plot_cells(cds, color_cells_by = "celltype",
           show_trajectory_graph = F,
           cell_size = 1)

# show trajectory
plot_cells(cds, color_cells_by = "celltype",
           show_trajectory_graph = T, label_leaves = F, label_branch_points = F, label_roots = F,
           cell_size = 1, group_label_size = 4)
plot_cells(cds, color_cells_by = "cluster",
           show_trajectory_graph = T, label_leaves = F, label_branch_points = F, label_roots = F,
           cell_size = 1, group_label_size = 4)

# plot expression of some random genes
plot_cells(cds, genes=c("Sox17"
                        ,"Hand2"),
           cell_size = 1)

plot_cells(cds, 
           show_trajectory_graph = T,color_cells_by = 'pseudotime', 
           label_leaves = F, label_branch_points = F, label_roots = F,label_groups_by_cluster =T,
           trajectory_graph_color = 'black', cell_size = 1)
plot_cells(cds, 
           show_trajectory_graph = T,color_cells_by = 'day', 
           label_leaves = F, label_branch_points = F, label_roots = F,label_groups_by_cluster =T,
           trajectory_graph_color = 'black', cell_size = 1)
plot_cells(
  cds = cds,
  color_cells_by = "pseudotime",
  show_trajectory_graph = TRUE, cell_size = 2
)
