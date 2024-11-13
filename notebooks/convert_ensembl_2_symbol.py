import scanpy as sc
import decoupler as dc

# Only needed for processing
import numpy as np
import pandas as pd
from anndata import AnnData


annot = sc.queries.biomart_annotations("hsapiens",
        ["ensembl_gene_id", "external_gene_name"],
        use_cache=False
    ).set_index("ensembl_gene_id")


def filter_nameless_genes(ensembl_names):
    return [item for item in ensembl_names if item.split('.')[0] in annot.index]
def convert_to_symbol(ensembl_names):
    ensembl_names = [item.split('.')[0] for item in ensembl_names]
    symbol_names = [annot.loc[ensembl_id,'external_gene_name'] for ensembl_id in ensembl_names]

    return symbol_names

def convert_file_to_symbol():
    pass