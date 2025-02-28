import pandas as pd
import numpy as np
import os
## print pwd
os.system('pwd')
import sys
sys.path.append('/home/amore/work')
from mirkitten import get_DE_dds_dict
import yaml



class DDS:
    """
    Class to prepare the DDS input of mirKat.

    It will take the individual DDSs of the different comparisons and will combine them into a single DataFrame.
    
    
    """
    def __init__(self, dds_files_path='../data/dds_files.yml', pvalue=0.05, threshold=0, interest='stat'):
        """
        Initializes the DDS class by loading the DDS files specified in a YAML file.
        """
        self.dds_dict = self._load_dds_files(dds_files_path)
        self.pvalue = pvalue
        self.threshold = threshold
        self.interest = interest
        self.dds_df = self.combine_dds_interest_genes(interest=interest)
    
    def _load_dds_files(self, dds_files_path):
        """
        Loads DDS files from a YAML file and returns a dictionary of DataFrames.
        """
        with open(dds_files_path, 'r') as file:
            dds_dict = yaml.safe_load(file)
        return {dds: pd.read_csv(file, index_col=0) for dds, file in dds_dict.items()}

    def combine_dds_interest_genes(self, interest='stat'):
        """
        Combines statistical values of interest from DDS comparisons into a single DataFrame.
        """    
        interest_dict = {}
        for comparison, dds_sig in self.dds_dict.items():
            stat = dds_sig[interest]
            interest_dict[comparison] = stat
        all_stat_df=None
        for comparison, stat in interest_dict.items():
            #comparison = comparison[:-5]
            if all_stat_df is None:
                all_stat_df = stat.rename(comparison).to_frame()
            else:
                all_stat_df = all_stat_df.join(stat.rename(comparison))
        all_stat_df.fillna(0, inplace=True)
        return all_stat_df
    
    def get_dds_df(self):
        """
        Returns the combined DDS DataFrame.
        """
        return self.dds_df
    def save_dds_df(self, path='results/dds_df.csv'):
        """
        Saves the combined DDS DataFrame to a CSV file.
        """
        self.dds_df.to_csv(path)
