import scanpy as sc
import decoupler as dc

# Only needed for processing
import numpy as np
import pandas as pd
from anndata import AnnData
from pydeseq2.dds import DeseqDataSet, DefaultInference
from pydeseq2.ds import DeseqStats


class DEA:
    """
    This class will be used to perform Differential Expression Analysis
    """
    def __init__(self, experiment='RNAseq_abundances_adjusted_combat_inmose_gene_symbol_expression.csv',
                  metadata_file='/home/amore/work/data/All_rna_samples_metadata.csv', 
                  save_folder= '/home/amore/work/results'):
        """
        
        """
        self.experiment = experiment
        self.metadata_file = metadata_file
        self.adata = pd.read_csv(f'/home/amore/work/data/{experiment}', index_col=0)
        self.metadata_df = pd.read_csv(metadata_file, index_col=0)
        self.inference = DefaultInference(n_cpus=8)
        if save_folder[-1] == '/':
            save_folder = save_folder[:-1]
        self.save_folder = save_folder
        self.sex = None

    def separate_by_sex(self, sex):
        """
        This function will separate
        """

        self.adata['Sex'] = self.metadata_df['Sex']
        self.adata = self.adata[self.adata['Sex']==sex]
        self.adata = self.adata.drop(columns='Sex')
        self.sex = sex

    def get_DEA(self, younger_group='young', older_group ='old', save=True):
        dds = self.get_dds()
        comparison = f'{younger_group}.vs.{older_group}'
        stat_res = DeseqStats(
            dds,
            contrast=["condition", younger_group, older_group],
            inference=self.inference
        )
        stat_res.summary()
        results_df = stat_res.results_df
        if save:
            if self.sex is not None:
                results_df.to_csv(f'{self.save_folder}/{self.experiment}_{comparison}_{self.sex}_DDS.csv', header=True)
            else:
                results_df.to_csv(f'{self.save_folder}/{self.experiment}_{comparison}_DDS.csv', header=True)
        return results_df
    @staticmethod
    def get_condition_mapping(ages, young=35, old=65):
        """
        This function will return a dictionary with the condition mapping
        """
        ages = ages.apply(lambda age: 'young' if age < young else ('old' if age >= old else 'middle'))
        return ages

    def get_dds(self):
        """
        This function will return the DESeq2 object dds needed to perfome the DEA
        """
        # Obtain genes that pass the thresholds
        ages=self.adata['Age']
        adata=self.adata.drop('Age',axis=1)
        conditions = self.get_condition_mapping(ages)
        andata=AnnData(adata.to_numpy(), obs=pd.DataFrame(ages))
        andata.obs_names = adata.index
        andata.var_names = adata.columns
        andata.obs['condition'] = conditions
        genes = dc.filter_by_expr(andata, group='condition', min_count=10, min_total_count=15, large_n=1, min_prop=1)
        andata = andata[:, genes].copy()

        # Build DESeq2 object
        
        dds = DeseqDataSet(
            adata=andata,
            design_factors='condition',
            refit_cooks=True,
            inference=self.inference,
        )
        dds.deseq2()

        return dds
    
    def save_DEA(self,
                 results_df,
                 comparison):
        """
        This function will save the results of the DEA
        """
        results_df.to_csv(f'{self.save_folder}/{self.experiment}_{comparison}.csv', header=True)



#### __main__ with no parameters, everything runs as default value. Creates object DEA and runs get_DEA ####
if __name__ == '__main__':
    dea = DEA()
    df = dea.get_DEA()
    print(df.head())

