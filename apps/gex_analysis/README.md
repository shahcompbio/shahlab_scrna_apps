# GEX_ANALYSIS

The main objective of this app is to do basic QC and pre-processing on scRNA Seq samples.

# NOTE

## Sample Level

The final h5ad has the following additional fields:
- obs.aliquot_id - This is the aliquot id of the analysis. Descriptive of the site, but does not uniquely identify the target.
- obs.sample_id - This is meant to directly link back to the sample to which the analysis is running on.
- layers.counts - This is the original counts, without normalization and log applied

# Flow

## Sample Level

Steps
1) Read the cellranger 10x matrix
2) Identify the doublets using (scrublet)[https://github.com/swolock/scrublet] and assign the scores.
3) Identify QC Metrics for MT(mito) Genes
4) Generate the violin plots that show the cells that will be excluded in QC.
5) Remove doublets and preform QC based on metrics in the options.
6) Save raw counts in layer counts.
7) Normalize and log.
8) Get highly expressing top_n_genes, store in the anndata, but do not filter.
9) Preform analysis (PCA, Neighbors, UMAP, Leiden Clustering, Rank Gene Groups by Leiden)
10) Plot the UMAP and save this final h5ad

Extra step: Differentially Expressed Genes (DEG)
11) Get the: gene, adj_pval, log_fc (log fold changes subtype (cluster)


# Options

    --total_counts : default=1000, type=int, help="total countes per cell minimum"
    
    --scrublet_score : default=0.35, type=click.FloatRange(min=0.0, max=1.0 help="scrublet: max threshold for the doublet score between [0, 1]"
    
    --n_prin_comps : default=10, type=int, help="scrublet: number of principal components used to embed the transcriptomes prior to k-nearest-neighbor graph construction"
    
    --pct_mito : default=25, type=int, help="the percentage of counts in mitochondrial genes maximum"
    
    --min_genes : default=500, type=int, help="minimum number of genes expressed required for a cell to pass filtering"
    
    --min_cells : default=3, type=int, help='minimum number of cells expressed required for a gene to pass filtering'
    
    --n_top_genes : default=5000, type=int, help='number of highly variable genes to filter for'

    -> --leiden_resolution : default=1.0, type=float, help='controls the coarseness of clustering'
    -> --subsig_pthresh : default=0.001, type=float, help='deg: the p threshold to check that the pval less than'
    -> --subsig_fctresh : default=0.5, type=float, help='deg: check that logfoldchanges is greater than fctresh'
    -> --subsig_top_n : default=50, type=int, help='deg: how many significant genes in subset to return'
    -> --subsig_sort_by : default=2, type=int, help='deg: sort based on negative of: 0-> gene, 1-> pval, 2->logfoldchanges

