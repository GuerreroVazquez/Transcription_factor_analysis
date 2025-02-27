import scanpy as sc
import decoupler as dc
import os
# Only needed for processing
import numpy as np
import pandas as pd
from anndata import AnnData



def convert_from_columns(experiment, species='hsapiens'):
    adata = pd.read_csv(f'/home/amore/work/data/{experiment}', columns=0)
 
    ensembl_names = adata.columns
    # remove ensembl_names that do not start with ENSG
    ensembl_names = [name for name in ensembl_names if name.startswith('ENSG')]
    ensembl_names = ensembl_names.to_list()
    symbols = convert_ensembl_symbol(ensembl_names, species)
    # For each column name that starts in ENSG, asign the corresponding gene symbol
    original_columns = adata.columns
    for column in original_columns:
        if column in ensembl_names:
            adata.rename(columns={column: symbols[ensembl_names.index(column)]}, inplace=True)
    return adata

def convert_from_rows(experiment, species='hsapiens'):
    adata = pd.read_csv(f'/home/amore/work/data/{experiment}', index_col=0)
    ensembl_names = adata.index
    # remove ensembl_names that do not start with ENSG
    ensembl_names = [name for name in ensembl_names if name.startswith('ENSG')]
    ensembl_names = ensembl_names.to_list()
    symbols = convert_ensembl_symbol(ensembl_names, species)
    # For each row name that starts in ENSG, asign the corresponding gene symbol
    original_rows = adata.index
    for row in original_rows:
        if row in ensembl_names:
            adata.rename(index={row: symbols[ensembl_names.index(row)]}, inplace=True)
    return adata


def convert_ensembl_symbol(ensembl_names, species='hsapiens'):
    if species=='hsapiens' and "anot.pkl" in os.listdir('/home/amore/work/data'):
        annot = pd.read_pickle('/home/amore/work/data/anot.pkl')
    else:
        annot = sc.queries.biomart_annotations(species,
                                               ["ensembl_gene_id", 
                                                "external_gene_name"],
                                                use_cache=False.set_index("ensembl_gene_id"))
    ensembl_names = [item for item in ensembl_names if item.split('.')[0] in annot.index]
    ensembl_names = [element.split('.')[0] for element in ensembl_names]
    symbols = [annot.loc[ensembl_id,'external_gene_name'] for ensembl_id in ensembl_names]
    return symbols






