import scrublet
import scanpy as sc
import pandas as pd
import numpy as np
import argparse
import json
from matplotlib import pyplot as plt
import matplotlib as mpl
from seaborn import violinplot, stripplot
import click

from os.path import join, dirname
from shutil import copyfile, SameFileError

plt.tight_layout()
mpl.rcParams['figure.dpi'] = 300

def identify_doublets(adata, n_prin_comps, scrublet_score):
    scrub = scrublet.Scrublet(adata.X)
    doublet_scores, _ = scrub.scrub_doublets(n_prin_comps=n_prin_comps)
    adata.obs["doublet_score"] = doublet_scores
    adata.obs["likely_singlet"] = adata.obs["doublet_score"].apply( 
            lambda score: "Singlet" if score < scrublet_score else "Doublet"
        )

def save_qc_graphs(adata, violin_plots_pre_filtering, min_genes, min_cells, total_counts, pct_mito ):
    fig, axs = plt.subplots(2, 2, figsize=(8,8))
    axs = axs.flatten()
    sc.pl.violin(
        adata,
        keys="n_genes_by_counts",
        groupby='likely_singlet',
        ax=axs[0],
        log=True
    )
    sc.pl.violin(
        adata,
        keys="total_counts",
        groupby='likely_singlet',
        ax=axs[1],
        log=True
    )
    sc.pl.violin(
        adata,
        keys="pct_counts_mt",
        groupby='likely_singlet',
        ax=axs[2],
        log=True
    )

    # Plot the genes - scanpy only plots obs
    axs[3].set(yscale="log")
    violinplot(data=adata.var['n_cells_by_counts'], ax=axs[3])
    stripplot(data=adata.var['n_cells_by_counts'], ax=axs[3], color='black', edgecolor='gray', size=1)

    # Get max for cutoffs
    max_genes = np.max(adata.obs['n_genes_by_counts'][ adata.obs['likely_singlet'] == 'Doublet' ])
    max_counts = np.max(adata.obs['total_counts'][ adata.obs['likely_singlet'] == 'Doublet' ])
    max_counts_mt = np.max(adata.obs['pct_counts_mt'][ adata.obs['likely_singlet'] == 'Doublet' ])

    # Color red everything that we will not take

    axs[0].axhline(min_genes, color='red')
    axs[0].fill_between([-0.5,1.5],[min_genes],[0], facecolor='red', alpha=0.2)
    axs[0].axvline(0.5, color='red')
    axs[0].fill_between([-0.5,0.5],[max_genes],[0], facecolor='red', alpha=0.2)
    axs[0].set_ylabel('')
    axs[0].set_title('n_genes_by_counts')

    axs[1].axhline(total_counts, color='red')
    axs[1].fill_between([-0.5,1.5],[total_counts],[0], facecolor='red', alpha=0.2)
    axs[1].axvline(0.5, color='red')
    axs[1].fill_between([-0.5,0.5],[max_counts],[0], facecolor='red', alpha=0.2)
    axs[1].set_ylabel('')
    axs[1].set_title('total_counts')

    axs[2].axhline(pct_mito, color='red')
    axs[2].fill_between([-0.5,1.5],[pct_mito],[100], facecolor='red', alpha=0.2)
    axs[2].axvline(0.5, color='red')
    axs[2].fill_between([-0.5,0.5],[max_counts_mt],[0], facecolor='red', alpha=0.2)
    axs[2].set_ylabel('')
    axs[2].set_title('pct_counts_mt')

    axs[3].axhline(min_cells, color='red')
    axs[3].fill_between([-0.5,0.5],[min_cells],[0], facecolor='red', alpha=0.2)
    axs[3].set_ylabel('')
    axs[3].set_title('genes_num_cells')

    fig.savefig(violin_plots_pre_filtering)


def process_deg(adata, output_deg, subsig_pthresh, subsig_fctresh, subsig_top_n, subsig_sort_by):
    # Split the new rank_genes data into easier to work with format
    genes = list(zip(*adata.uns['rank_genes_groups']["names"].tolist()))
    subtypes = adata.uns['rank_genes_groups']["names"]
    subtypes = subtypes.dtype.names
    adjpvals = list(zip(*adata.uns['rank_genes_groups']["pvals_adj"].tolist()))
    logfc = list(zip(*adata.uns['rank_genes_groups']["logfoldchanges"].tolist()))

    def subset_significant(genes, pvalues, fcs, 
                            pthresh=subsig_pthresh,
                            fcthresh=subsig_fctresh,
                            top_n=subsig_top_n,
                            sort_by=subsig_sort_by
                        ):
        # Do the actual threshold test and return only the genes of the subset that pass
        subset = [
            (gene, pvalue, fc)
            for gene, pvalue, fc in zip(genes, pvalues, fcs)
            if pvalue < pthresh and fc > fcthresh
        ]

        subset.sort(key=lambda x:-1.0 *x[sort_by])
        return subset[:top_n]

    df_genes   = []
    df_pvals   = []
    df_fcs     = []
    df_subtype = []

    for subtype, gene, pval, fc in zip(subtypes, genes, adjpvals, logfc):
        subset = subset_significant(gene, pval, fc)
        df_genes += [x[0] for x in subset]
        df_pvals += [x[1] for x in subset]
        df_fcs += [x[2] for x in subset]
        df_subtype += [subtype for _ in subset]

    data = {"gene":df_genes, "adj_pval":df_pvals, "log_fc": df_fcs, "subtype":df_subtype}
    df = pd.DataFrame.from_dict(data)

    df.to_csv(output_deg, sep="\t", index=False)

def run_gex_analysis(params: dict):
    """
    Params file specific in `run`
    """
    # Process Mandatory Args
    system_id = params['system_id']
    aliquot_id = params['aliquot_id']
    cellranger_matrix = params['cellranger_matrix'] # Folder with matrix/barcode/features .tsv.gz
    output_dir = params['output_dir']

    # Process Optional Args (All should have default values)
    total_counts = params['total_counts']
    scrublet_score = params['scrublet_score']
    pct_mito = params['pct_mito']
    min_genes = params['min_genes']
    min_cells = params['min_cells']
    n_prin_comps = params['n_prin_comps']
    leiden_resolution = params['leiden_resolution']
    subsig_pthresh = params['subsig_pthresh']
    subsig_fctresh = params['subsig_fctresh']
    subsig_top_n = params['subsig_top_n']
    subsig_sort_by = params['subsig_sort_by']
    n_top_genes = params['n_top_genes']

    # Outputs to save
    output_h5ad = join(output_dir, "output_h5ad.h5ad")
    output_deg = join(output_dir, "output_deg.tsv")
    umap_url = join(output_dir, "umap.png")
    violin_plots_pre_filtering = join(output_dir, "violin_plots_pre_filtering.png")
    highly_variable_image = (output_dir, "highly_variable_image.png")

    # If we are to load based on CELLRANGER data
    #Load 10X Gex
    adata = sc.read_10x_mtx(cellranger_matrix)

    #Identify Doublet Cells
    identify_doublets(adata, n_prin_comps, scrublet_score)

    # Get mitocondrial
    adata.var['mt'] = adata.var_names.str.startswith('MT-')
    sc.pp.calculate_qc_metrics(adata, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)

    # Get some QC Graphs before we do any filtering so that we know where the cutoffs are
    save_qc_graphs(adata, violin_plots_pre_filtering, min_genes, min_cells, total_counts, pct_mito )

    # Remove the doublets
    adata = adata[adata.obs["doublet_score"] < scrublet_score]

    #QC filters
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)
    adata = adata[adata.obs.total_counts > total_counts, ]
    adata = adata[adata.obs.pct_counts_mt < pct_mito, ]

    # Normalize and log, but save original counts in layer
    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata)
    sc.pp.log1p(adata)

    # Get the highest expressing genes
    sc.pp.highly_variable_genes(
        adata,
        n_top_genes=n_top_genes,
    )

    # Save figure of highly expressing genes BUT DO NOT CUT OFF
    sc.pl.highly_variable_genes(
        adata,
        log=True,
        show=False,
        save=highly_variable_image
    )

    #Generate UMAP embedding
    sc.pp.pca(adata, svd_solver="arpack")
    sc.pp.neighbors(adata)
    sc.tl.umap(adata)

    # Generate Leiden clusters and rank genes
    sc.tl.leiden(adata, resolution = leiden_resolution)
    sc.tl.rank_genes_groups(adata, 'leiden', method='wilcoxon')

    # Prefix Leiden rows with Cluster_
    adata.obs['leiden'] = 'Cluster_' + adata.obs['leiden'].astype(str)

    # Now get the umap data
    umap1, umap2 = list(zip(*adata.obsm['X_umap'])) # UMAP 1 and UMAP 2
    adata.obs = adata.obs.assign(UMAP_1=umap1, UMAP_2=umap2 )

    # Plot the umap
    fig = plt.figure()
    ax = fig.subplots(1, 1)
    sc.pl.umap(adata, color=['leiden'], use_raw=False, ax=ax)
    plt.tight_layout()
    fig.savefig(umap_url)

    # Save the aliquot_id and sample_id to the h5ad
    adata.obs['system_id'] = system_id
    adata.obs['aliquot_id'] = aliquot_id

    #### H5AD WRITE HERE
    # Write out our currnt h5da
    adata.write(output_h5ad)

    #### DEG PROCESSING HERE
    process_deg(adata, output_deg, subsig_pthresh, subsig_fctresh, subsig_top_n, subsig_sort_by)

@click.group(invoke_without_command=True)
@click.pass_context
def cli(ctx):
    if ctx.invoked_subcommand is None:
        click.echo(
            "For help, please run: `gex_analysis --help`"
        )

@cli.command()
@click.argument("save_path")
def get_default_config(save_path):
    """
    Get the default config to be specific as params file for run.
    - save_path
        Where to save the default json configuration
    """
    try:
        default_config = join(dirname(__file__), 'default_config.json')
        copyfile(default_config, save_path)
    except SameFileError as SFE:
        click.echo("Can't copy same file! If the path the same as the installed app path?")
    except Exception as e:
        raise e

@cli.command()
@click.argument("params_file")
def run(params_file):
    """
    Run the gex_analysis app.

    Parameters:
    - params_file
        Configuration JSON with the properties listed below.

    The params file has to have the following options in json:
    - cellranger_matrix
        The output path from `cellranger count`
    - output_dir
        Where to save the analyses
    - system_id
        The exact ID to use to reference the sample (can be empty)
    - aliquot_id
        A more descriptive ID to use to explain the experiment (can be empty)

    Optional parameters to consider:
    - total_counts
        The minimum total counts for each cell
    - scrublet_score
        The min scrublet score to be considered a singlet
    - pct_mito
        The maximum percent mito genes to keep
    - min_genes
        The minimum number of genes expressed in a cell to be kept
    - min_cells
        The minimum number of cells per gene to be kept
    - n_prin_comps
        For scrublet - Number of principal components to use for embedding cells prior to
        building kNN graph.
    - leiden_resolution
        Leiden resolution when doing leiden clustering
    - subsig_pthresh
    - subsig_fctresh
    - subsig_top_n
    - subsig_sort_by
    - n_top_genes

    """
    # Read JSON file
    with open(params_file) as opt_file:
        params = json.load(opt_file)

    run_gex_analysis(params)

if __name__ == '__main__':
    cli()
