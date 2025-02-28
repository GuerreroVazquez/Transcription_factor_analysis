import os
import pandas as pd
import decoupler as dc
## print pwd
os.system('pwd')
import sys
sys.path.append('/home/amore/work')
from mirkitten import get_DE_genes_df
import yaml



class Pathways:
    def __init__(self, dds_dict, sel_db = ['go_molecular_function',
                                  'go_cellular_component',
                                  'go_biological_process',
                                  'reactome_pathways',
                                  'kegg_pathways', 'hallmark'], pvalue=0.05, threshold=None, interest='stat', pathway_pvalue=0.05):
        if pathway_pvalue is None:
            pathway_pvalue=0.05
        if "msigdb" in os.listdir('/home/amore/work/data'):
            msigdb = pd.read_csv('msigdb.csv')
        else:
            msigdb = dc.get_resource('MSigDB')
        self.msigdb =msigdb[msigdb['collection'].isin(sel_db)]
        self.msigdb = self.msigdb[~self.msigdb.duplicated(['geneset', 'genesymbol'])]
        self.dds_dict = self._load_dds_files(dds_files_path=dds_dict)
        self.sel_db = sel_db
        self.pvalue = pvalue
        self.threshold = threshold
        self.interest = interest
        self.pathway_pvalue = pathway_pvalue
        self.enriched_ORA_pathway = None
    def _load_dds_files(self, dds_files_path):
        """
        Loads DDS files from a YAML file and returns a dictionary of DataFrames.
        """
        with open(dds_files_path, 'r') as file:
            dds_dict = yaml.safe_load(file)
        return {dds: pd.read_csv(file, index_col=0) for dds, file in dds_dict.items()}
    def set_sel_db(self, sel_db):
        self.sel_db = sel_db
    def load_sel_df_from_file(self, path):
        """
        This function will read a text file that has the selected databases
        and take them a list of string that then will be added to self.sel_db
        """
        with open(path, 'r') as file:
            self.sel_db = file.readlines()
    
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
        self.msigdb = self.msigdb[~self.msigdb.duplicated(['geneset', 'genesymbol'])]
        degs = [item for item in degs if "Unnamed" not in item and "NAN" not in item]
        self.enr_pvals = dc.get_ora_df(
            df=degs,
            net=self.msigdb,
            source='geneset',
            target='genesymbol'
        )
        print(self.enr_pvals.head())
        return self.enr_pvals

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
            print(pathways.head())
            pathways_dict[comparison] = pathways

        all_scores_pathway_df=None
        for comparison, score in pathways_dict.items():
            if all_scores_pathway_df is None:
                all_scores_pathway_df = score.rename(comparison).to_frame()
            else:
                all_scores_pathway_df = all_scores_pathway_df.join(score.rename(comparison))
        all_scores_pathway_df.fillna(0, inplace=True)
        return all_scores_pathway_df
    def save_enriched_pathways_combined_ORA(self, path='results/enriched_pathways_combined_ORA.csv'):
        pathway_df = self.get_enriched_pathways_combined_ORA()
        pathway_df.to_csv(path)
    def set_enriched_ORA_pathway(self ):
        self.enriched_ORA_pathway = self.get_enriched_pathways_combined_ORA()
        #print (self.enriched_ORA_pathway.head())
        
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
            if pathway in msigdb.index:
                genes = msigdb.loc[pathway, 'genesymbol']
                gene_list.append(genes)
            else: 
                pass#print (pathway)
        enrriched_pathway_df['genes'] = gene_list
        
        return enrriched_pathway_df

    def get_enriched_pathways_with_genes_ORA(self):
        """
        This function will return the enriched pathways with the genes
        """
        if self.enriched_ORA_pathway is None:
            self.set_enriched_ORA_pathway()
        return self.attach_gene_list_to_pathway()
    
    def save_enrriched_pathways_with_genes_ORA(self, path):
        df = self.get_enriched_pathways_with_genes_ORA()
        df.to_csv(path)


