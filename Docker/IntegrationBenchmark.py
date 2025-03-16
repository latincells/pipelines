import numpy as np
import scanpy as sc
from scib_metrics.benchmark import Benchmarker ############
import pandas as pd
import datashader as ds ############  install conda-forge::datashader
import datashader.transfer_functions as tf ############
import datashader.bundling as bd ############
import matplotlib.pyplot as plt
import colorcet ############ install conda-forge::colorcet
import matplotlib.colors
import matplotlib.cm

import bokeh.plotting as bpl ############ install conda-forge::bokeh
import bokeh.transform as btr ############
import holoviews as hv ############ install conda-forge::holoviews
import holoviews.operation.datashader as hd ############

#################################################################################################
import scipy  # Import the scipy library, which provides functions for scientific computing, including handling sparse matrices.
import scanorama  # Import the Scanorama library for batch effect correction and integration of single-cell RNA-seq data.
import numpy as np  # Import the numpy library for numerical operations, particularly for handling arrays.

from harmony import harmonize  # Import the 'harmonize' function from the Harmony library, which is used for correcting batch effects in PCA-reduced data.
from rich import print

import scvi
#################################################################################################
Dir = "LatinCells_Results/Chile_Tests/4-Seurat_QC-Integration/DataForCelltypist/"
adata = sc.read_10x_mtx(Dir)
metadata = pd.read_csv( Dir + 'Metadata.csv',sep="\t",low_memory=False)
adata.obs = adata.obs.join(metadata)
#umap = pd.read_csv( Dir + 'Cell.embeddings.UMAP.csv',sep="\\t",low_memory=False)
#pca = pd.read_csv( Dir + 'Cell.embeddings.PCA.csv',sep="\\t",low_memory=False)
#################################################################################################
# 4. Exploratory Data Analysis
adata.obs["orig.ident"].value_counts()
adata.obs["SNG.BEST.GUESS"].value_counts()
adata.obs["Country"].value_counts()
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
sc.pp.highly_variable_genes(adata, n_top_genes=2000, flavor="cell_ranger", batch_key="orig.ident")
sc.tl.pca(adata, n_comps=30, use_highly_variable=True)

print(adata.shape)
adata = adata[:, adata.var.highly_variable].copy()
print(adata.shape)

#################################################################################################
# 5. Data Visualization

# Performs Principal Component Analysis (PCA) on the selected genes
sc.tl.pca(adata)

# Saves the PCA results in the .obsm attribute under the key "Unintegrated"
# This essentially creates a copy of the PCA results, labeling it as "Unintegrated"
# to distinguish it from other potential embeddings or analyses.
adata.obsm["Unintegrated"] = adata.obsm["X_pca"]

# Computes the nearest neighbors graph for the dataset.
# This is a crucial step for many downstream analyses, such as UMAP or clustering,
# as it defines the local neighborhood of each cell based on the previously computed PCA components.
sc.pp.neighbors(adata, use_rep = "Unintegrated")

# Applies UMAP (Uniform Manifold Approximation and Projection) to the dataset.
# UMAP is a dimensionality reduction technique that projects the high-dimensional data
# into a 2D or 3D space, making it easier to visualize clusters and patterns among the cells.
sc.tl.umap(adata)

# Generates a UMAP plot to visualize the data with specific coloring based on cell type.
sc.pl.umap(
    adata,
    color=["orig.ident"],  # Colors the cells in the UMAP plot according to the "cell_type" metadata.
                          # You can change "cell_type" to other variables like "batch" or "patientGroup"
                          # to observe the effects of integration or to highlight different metadata categories.
    # legend_loc='on data', # Optional: Uncomment to place the legend directly on the plot.
    title="Unintegrated dataset",  # Sets the title of the plot to "Unintegrated dataset".
    frameon=False,  # Removes the frame around the plot for a cleaner appearance.
    ncols=1,  # Specifies that the plot should be displayed in a single column layout.
    save = "_Unintegrated_dataset.orig.ident.png"
)

sc.pl.umap(
    adata,
    color=["SNG.BEST.GUESS"],
    title="Unintegrated dataset",
    frameon=False,
    ncols=1,
    save = "_Unintegrated_dataset.BEST.GUESS.png"
)

#################################################################################################
# 6. Integration and Benchmarking
# Integration with Scanorama:

# List of adata per batch
batch_cats = adata.obs["orig.ident"].cat.categories  # Get the unique categories (batch identifiers) from the 'batch' column of the 'adata' object.
adata_list = [adata[adata.obs["orig.ident"] == b].copy() for b in batch_cats]  # Create a list of 'adata' objects, one for each batch, by filtering the original 'adata' based on the batch identifier.

# Convert to CSR format if needed
for i, ad in enumerate(adata_list):  # Iterate over each 'adata' object in the 'adata_list'.
    if isinstance(ad.X, scipy.sparse.csc_matrix):  # Check if the data matrix 'X' is in CSC (Compressed Sparse Column) format.
        ad.X = ad.X.tocsr()  # Convert the matrix from CSC to CSR (Compressed Sparse Row) format if necessary. CSR is often preferred for certain operations.

# Integrate with Scanorama
scanorama.integrate_scanpy(adata_list)  # Use Scanorama to integrate the list of 'adata' objects to correct for batch effects and combine the data into a unified representation.

# Prepare the integrated matrix
adata.obsm["Scanorama"] = np.zeros((adata.shape[0], adata_list[0].obsm["X_scanorama"].shape[1]))  # Initialize a new array in 'adata.obsm' to store the integrated matrix from Scanorama. The array has the same number of rows as the original 'adata' and columns equal to the Scanorama matrix dimension.
for i, b in enumerate(batch_cats):  # Iterate over each batch category.
    adata.obsm["Scanorama"][adata.obs["orig.ident"] == b] = adata_list[i].obsm["X_scanorama"]  # Assign the integrated Scanorama matrix for each batch to the corresponding entries in the new matrix.

# harmony
# Perform batch effect correction using Harmony
adata.obsm["Harmony"] = harmonize(adata.obsm["X_pca"], adata.obs, batch_key="orig.ident")

# Integration with scVI
adata.layers['counts'] = adata.X
# Store the raw count matrix (adata.X) in the 'counts' layer of the AnnData object.
# This allows the SCVI model to use this data as the input for its analysis.

scvi.model.SCVI.setup_anndata(adata, layer="counts", batch_key="orig.ident")
# Prepare the AnnData object for SCVI by specifying that the raw counts are stored in the 'counts' layer.
# The 'batch' column in 'adata.obs' is used to account for batch effects during model training.

vae = scvi.model.SCVI(adata, gene_likelihood="nb", n_layers=2, n_latent=30)
# Initialize the SCVI model with the following parameters:
# - 'gene_likelihood="nb"' specifies that the negative binomial distribution is used for modeling gene expression.
# - 'n_layers=2' sets the number of hidden layers in the neural network.
# - 'n_latent=30' sets the number of latent dimensions to learn from the data.

vae.train()
# Train the SCVI model on the AnnData object. This step fits the model to the data and learns the latent representations.

adata.obsm["scVI"] = vae.get_latent_representation()
# Store the learned latent representations from the SCVI model into 'adata.obsm' under the key 'scVI'.
# This representation captures the underlying structure of the data after accounting for batch effects and other variations.

# Integration with scANVI
lvae = scvi.model.SCANVI.from_scvi_model(
     vae,
     adata=adata,
     labels_key="cell_type",
     unlabeled_category="Unknown",
 )
# Initialize the SCANVI model using the previously trained SCVI model.
# 'vae' is the trained SCVI model that serves as the starting point.
# 'adata' is the AnnData object containing the data.
# 'labels_key="cell_type"' specifies that the cell type labels are stored in the 'cell_type' column of 'adata.obs'.
# 'unlabeled_category="Unknown"' designates the category for unlabeled cells in the dataset.

lvae.train(max_epochs=20, n_samples_per_label=100)
# Train the SCANVI model with the specified parameters:
# - 'max_epochs=20' sets the maximum number of training epochs.
# - 'n_samples_per_label=100' specifies the number of samples to use per label during training.

adata.obsm["scANVI"] = lvae.get_latent_representation()
# Store the learned latent representations from the SCANVI model into 'adata.obsm' under the key 'scANVI'.
# This representation captures cell-type-specific latent features after incorporating the cell type labels and learning from both labeled and unlabeled cells.

# Dimensionality Reduction with Integration Neighbors

sc.pp.neighbors(adata, use_rep="Scanorama")
sc.tl.umap(adata)
sc.pl.umap(
    adata,
    color=["orig.ident"],
    title="Scanorama",
    frameon=False,
    ncols=1,
    save = "_Scanorama.png"
)

sc.pp.neighbors(adata, use_rep="Harmony")
sc.tl.umap(adata)
sc.pl.umap(
    adata,
    color=["orig.ident"],
    title="Harmony",
    frameon=False,
    ncols=1,
    save = "_Harmony.png"
)

sc.pp.neighbors(adata, use_rep="scVI")
sc.tl.umap(adata)
sc.pl.umap(
    adata,
    color=["orig.ident"],
    title="scVI",
    frameon=False,
    ncols=1,
    save = "_scVI.png"
)

# Initialize the Benchmarker object to evaluate the performance of different embeddings.
bm = Benchmarker(
    adata,  # The AnnData object containing the single-cell data and various embeddings.
    batch_key="orig.ident",  # The key in 'adata.obs' that holds batch information for evaluating batch effects.
    label_key="cell_type",  # The key in 'adata.obs' that contains cell type labels for assessing clustering performance.
    embedding_obsm_keys=["Unintegrated", "Scanorama", "Harmony", "scVI"], # , "scANVI" # The keys in 'adata.obsm' corresponding to the embeddings to be compared.
    n_jobs=2,  # The number of parallel jobs to use for computation, which can help speed up the benchmarking process.
)

# Run the benchmarking process to compare the performance of the specified embeddings.
# This involves evaluating how well each embedding method handles batch effects and clustering accuracy, based on the provided batch and label information.
bm.benchmark()
bm.plot_results_table()

from rich import print

df = bm.get_results(min_max_scale=False)
print(df)
df.transpose()

end_time = time.time()
elapsed_time = end_time - start_time

print(f'Total time taken: {elapsed_time:.2f} seconds')