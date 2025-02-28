import pandas as pd
import numpy as np
from kneed import KneeLocator

def get_DE_genes_df(dds:pd.DataFrame, pvalue:float=0.05, threshold=None, interest='stat'):
    """
    This function will take a dataframe with the dds results and will return a dataframe with the DE genes
    """
    if threshold is None:
        threshold = get_elbow_point_threshold(df = dds, interest=interest)
    return dds[(dds['padj'] < pvalue) & (dds[interest] > threshold)]

def get_elbow_point_threshold(df:pd.DataFrame, interest='stat'):
    """
    This function will return the threshold for the elbow point
    """

    # Step 1: Sort by evaluation (ascending or descending)
    df = df.sort_values(by=interest, ascending=False).reset_index(drop=True)

    # Step 2: Normalize values (optional)
    y = np.arange(len(df))  # Indices (gene positions)
    x = df[interest].values

    # Step 3: Apply Kneedle algorithm to find the elbow point
    kneedle = KneeLocator(x, y, curve='convex', direction='decreasing')  # Direction depends on sorting
    elbow_point = kneedle.elbow
    return elbow_point


def get_DE_dds_dict(dds_dict:dict, pvalue:float=0.05, threshold=None, interest='stat'):
    """
    This function will take a dictionary with the dds results and will return a dictionary with the DE genes
    """
    de_dict = {}
    for comparison, dds in dds_dict.items():
        de_dict[comparison] = get_DE_genes(dds, pvalue=pvalue, threshold=threshold, interest=interest)
    return de_dict

def get_DE_genes(dds:pd.DataFrame, pvalue:float=0.05, threshold=None, interest='stat'):
    """
    This function will take a dataframe with the dds results and will return a dataframe with the DE genes
    """
    de_df = get_DE_genes_df(dds, pvalue=pvalue, threshold=threshold, interest=interest)
    return list(de_df.index)