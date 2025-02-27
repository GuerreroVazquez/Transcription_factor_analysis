import pandas as pd
import yaml
from kneed import KneeLocator
import numpy as np
import os
import decoupler as dc

def get_DDS(dds_files_path = '../data/dds_files.yml'):

    """
    This funtion will take a yaml file, the yml file should have the following structure:
    dds1: /path/to/dds1.csv
    dds2: /path/to/dds2.csv
    etc
    The function will read the csv files and return a dictionary with the dataframes of the dds files.
    Important to mention that each csv file itself should have the columns
     baseMean  log2FoldChange     lfcSE      stat    pvalue      padj

    And the index should be the gene names in gene_symbol format. If there is not, please use the convert_ensembl2symbol function

    :param dds_files_path: str, path to the yml file
    :return: dict, dictionary with the dataframes of the dds files. The keys are the dds names (or comparisons) 
        and the values are the dataframes od the DDS dataframe

    """
    with open(dds_files_path, 'r') as file:
        dds_dict = yaml.safe_load(file)
    df_dds_dict = {}
    for dds, file in dds_dict.items():
        df_dds_dict[dds] =  pd.read_csv(file, index_col=0)
    return df_dds_dict


def get_DE_genes_df(dds:pd.DataFrame, pvalue:float=0.05, threshold=None, interest='stat'):
    """
    This function will take a dataframe with the dds results and will return a dataframe with the DE genes
    """
    if threshold is None:
        threshold = get_elbow_point_threshold(df = dds, interest=interest)
    return dds[(dds['padj'] < pvalue) & (dds[interest] > threshold)]
def get_DE_genes(dds:pd.DataFrame, pvalue:float=0.05, threshold=None, interest='stat'):
    """
    This function will take a dataframe with the dds results and will return a dataframe with the DE genes
    """
    de_df = get_DE_genes_df(dds, pvalue=pvalue, threshold=threshold, interest=interest)
    return list(de_df.index)

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

def combine_dds_interest_genes(dds_dict:dict, interest='stat'):
    """
    This function will take a dictionary with the dds results and will return a dictionary with the genes
    """

    interest_dict = {}
    for comparison, dds_sig in dds_dict.items():
        stat = dds_sig[interest]
        interest_dict[comparison] = stat
    all_stat_df=None
    for comparison, stat in interest_dict.items():
        #comparison = comparison[:-5]
        if all_stat_df is None:
            all_stat_df = stat.rename(comparison).to_frame()
        else:
            all_stat_df = all_stat_df.join(stat.rename(comparison))
    return all_stat_df
    

class DDS:
    def __init__(self, dds_files_path='../data/dds_files.yml', pvalue=0.05, threshold=0, interest='stat'):
        """
        Initializes the DDS class by loading the DDS files specified in a YAML file.
        """
        self.dds_dict = self._load_dds_files(dds_files_path)
        self.pvalue = pvalue
        self.threshold = threshold
        self.interest = interest
    
    def _load_dds_files(self, dds_files_path):
        """
        Loads DDS files from a YAML file and returns a dictionary of DataFrames.
        """
        with open(dds_files_path, 'r') as file:
            dds_dict = yaml.safe_load(file)
        return {dds: pd.read_csv(file, index_col=0) for dds, file in dds_dict.items()}

    
    def get_DE_genes(self, dds: pd.DataFrame, pvalue: float = 0.05):
        """
        Returns a list of differentially expressed (DE) genes.
        """
        de_df = dds[(dds['padj'] < pvalue)]
        return list(de_df.index)

    def get_DE_dds_dict(self, pvalue: float = 0.05, threshold=None, interest='stat'):
        """
        Returns a dictionary of differentially expressed (DE) genes for all DDS comparisons.
        """
        return {comparison: self.get_DE_genes(dds, pvalue, threshold, interest) for comparison, dds in self.dds_dict.items()}

    def combine_dds_interest_genes(self, interest='stat'):
        """
        Combines statistical values of interest from DDS comparisons into a single DataFrame.
        """
        de_dict = get_DE_dds_dict(pvalue=self.pvalue, threshold=self.threshold, interest=self.interest)
        interest_dict = {comp: dds[interest] for comp, dds in de_dict.items()}
        all_stat_df = None
        for comp, stat in interest_dict.items():
            if all_stat_df is None:
                all_stat_df = stat.rename(comp).to_frame()
            else:
                all_stat_df = all_stat_df.join(stat.rename(comp))
        # fill NaN values with 0
        all_stat_df.fillna(0, inplace=True)
        return all_stat_df


class Enrichment:
    def __init__(self, dds_dict, sel_db = ['go_molecular_function',
                                  'go_cellular_component',
                                  'go_biological_process',
                                  'reactome_pathways',
                                  'kegg_pathways', 'hallmark'], pvalue=0.05, threshold=None, interest='stat', pathway_pvalue=0.05):
        if "msigdb" in os.listdir('/home/amore/work/data'):
            msigdb = pd.read_csv('msigdb.csv')
        else:
            msigdb = dc.get_resource('MSigDB')
        msigdb =msigdb[msigdb['collection'].isin(sel_db)]
        self.dds_dict = dds_dict
        self.sel_db = sel_db
        self.pvalue = pvalue
        self.threshold = threshold
        self.interest = interest
        self.pathway_pvalue = pathway_pvalue
    
    def select_only_specific_pathways(self,pathways_interest):
        """
        This function will take the pathways of interest
        """
        self.msigdb = self.msigdb[self.msigdb['geneset']==pathways_interest]
        pass
    def restart_msigdb(self):
        if "msigdb" in os.listdir('/home/amore/work/data'):
            self.msigdb = pd.read_csv('msigdb.csv')
        else:
            self.msigdb = dc.get_resource('MSigDB')
        self.msigdb =self.msigdb[self.msigdb['collection'].isin(self.sel_db)]
    def show_msigdb_genesets(self):
        return self.msigdb['geneset'].unique()
    def show_msigdb_collections(self):
        return self.msigdb['collection'].unique()
    
    def run_ora(self, degs):
        self.enr_pvals = dc.get_ora_df(
            df=degs,
            net=self.msigdb,
            source='geneset',
            target='genesymbol'
        )

    def get_enriched_pathways(self, df):
        """
        This function will takes the DE genes and will run the Pathway enrichment ORA.

        :param df: pd.DataFrame, dataframe with the DE genes
        :return: pd.DataFrame, dataframe with the enriched pathways
        """
        DE_df = get_DE_genes_df(df)
        degs = DE_df.index
        pathway_df = self.run_ora(degs)
        pathway_df = pathway_df[pathway_df['p-value'] < self.pathway_pvalue]
        enriched_pathways = pathway_df['Combined score']
        
        return enriched_pathways

    def get_enriched_pathways_combined_ORA(self):
        """
        This function will take the dictionary of the already DE genes, and will run the Pathway enrichment ORA.
        """
        pathways_dict = {}
        for comparison, dds in self.dds_dict.items():
            pathways = self.get_enriched_pathways(dds)
            pathways_dict[comparison] = pathways

        all_scores_pathway_df=None
        for comparison, score in pathways_dict.items():
            if all_scores_pathway_df is None:
                all_scores_pathway_df = score.rename(comparison).to_frame()
            else:
                all_scores_pathway_df = all_scores_pathway_df.join(score.rename(comparison))
        all_scores_pathway_df.fillna(0, inplace=True)
        return all_scores_pathway_df

    def set_enriched_ORA_pathway(self ):
        self.enriched_ORA_pathway = self.get_enriched_pathways_combined_ORA()
        
    def attach_gene_list_to_pathway(self, enrriched_pathway_df=None):
        """
        This function will take all the pathways in enriched_ORA_pathway,
          will look at the msigdb and will attach the gene list to the pathway in a new column "genes"
        """
        if enrriched_pathway_df is None:
            if self.enriched_ORA_pathway is None:
                self.set_enriched_ORA_pathway()
            enrriched_pathway_df = self.enriched_ORA_pathway

        # The msigdb has columns genesymbol, collection, geneset.
        # for each pathway in the enriched_ORA_pathway, we will look at the msigdb and filter by the geneset. 
        # It will take all the genes in genesymbol and will attach them to the pathway in a new column "genes" as a list
        msigdb = self.msigdb
        msigdb = msigdb.set_index('geneset')
        gene_list = []

        for pathway in enrriched_pathway_df.index:
            genes = msigdb.loc[pathway, 'genesymbol']
            gene_list.append(genes)
        enrriched_pathway_df['genes'] = gene_list
        
        return enrriched_pathway_df

    def get_enriched_pathways_with_genes(self):
        """
        This function will return the enriched pathways with the genes
        """
        if self.enriched_ORA_pathway is None:
            self.set_enriched_ORA_pathway()
        return self.attach_gene_list_to_pathway()
    


class TranFact():
    def __init__(self, dds_dict, pvalue=0.05, threshold=None, interest='stat'):
        self.dds_dict = dds_dict
        self.pvalue = pvalue
        self.threshold = threshold
        self.interest = interest
        self.colletri = self.set_colletri()
        self.mats = {}
    def set_colletri(self, species='human'):
        if species=='human' and "collectri" in os.listdir('/home/amore/work/data'):
            collectri = pd.read_csv('collectri.csv')
        else:
            collectri = dc.get_collectri(organism=species, split_complexes=True)
        self.collectri = collectri
        return collectri

    def get_mat(self, dds, interest='stat', comparison='comparison', experiment='experiment'):
        """
        This function will take a dataframe with the dds results and will return a dataframe with the DE genes
        """
        if f'mat_{experiment}_{comparison}.csv' in os.listdir('results'):
            mat = pd.read_csv(f'results/mat_{experiment}_{comparison}.csv', index_col=0)
            mat = mat.T
        else:
            mat = dds[[interest]].T.rename(index={interest: comparison})
            mat.T.to_csv(f'results/mat_{experiment}_{comparison}.csv')
        self.mat[comparison]= mat
        return mat
    def get_tf_acts(self, mat, experiment='experiment', comparison='comparison'):
        if f'tf_acts_{experiment}_{comparison}.csv' in os.listdir('results'):
            tf_df = pd.read_csv(f'results/tf_acts{experiment}_{comparison}.csv', index_col=0)
        else:
            tf_acts, tf_pvals = dc.run_ulm(mat=mat, net=self.collectri, verbose=True)
            tf_df = pd.DataFrame(tf_acts.T)
            tf_df['pvals']=tf_pvals.T
            tf_df.to_csv(f'results/tf_acts{experiment}_{comparison}.csv')

    def get_stuff(self, dds, interest='stat', comparison='comparison', experiment='experiment'):
        """
        This function will take a dataframe with the dds results and will return a dataframe with the DE genes
        """
        mat = self.set_mat(dds=dds, interest=interest, comparison=comparison, experiment=experiment)
        
        collectri = self.collectri

        
        #tf_df = pd.read_csv(f'results/tf_acts{experiment}_{comparison}.csv', index_col=0)
        tf_df = tf_df[tf_df['pvals']<0.1]
        tf_acts=pd.DataFrame(tf_df[comparison])
        tf_acts=tf_acts.T
        values = tf_acts.iloc[0]
        down_reg = values.sort_values(ascending=True)[:5].index.to_list()
        up_reg = values.sort_values(ascending=False)[:5].index.to_list()
        up_down_reg = down_reg.copy()
        up_down_reg.extend(up_reg)

