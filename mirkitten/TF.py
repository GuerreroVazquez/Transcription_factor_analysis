import os
import pandas as pd
import decoupler as dc
import yaml


class TranFact():
    def __init__(self, dds_files_path, pvalue=0.05, threshold=None, interest='stat', species='human'):
        dds_dict = self._load_dds_files(dds_files_path=dds_files_path)
        self.dds_dict = dds_dict
        self.pvalue = pvalue
        self.threshold = threshold
        self.interest = interest
        self.colletri = self.set_colletri(species=species)
        self.mats = {}
        self.tf_dfs = {}
    def _load_dds_files(self, dds_files_path):
        """
        Loads DDS files from a YAML file and returns a dictionary of DataFrames.
        """
        with open(dds_files_path, 'r') as file:
            dds_dict = yaml.safe_load(file)
        return {dds: pd.read_csv(file, index_col=0) for dds, file in dds_dict.items()}
    
    def set_colletri(self, species='human'):
        if species=='human' and "collectri" in os.listdir('/home/amore/work/data'):
            collectri = pd.read_csv('collectri.csv')
        else:
            collectri = dc.get_collectri(organism=species, split_complexes=True)
        self.collectri = collectri
        return collectri

    def get_mat(self, dds, interest='stat', comparison='comparison'):
        """
        This function will take a dataframe with the dds results and will return a dataframe with the DE genes
        """
        if comparison in self.mats:
            mat = self.mats[comparison]
        if f'mat_{comparison}.csv' in os.listdir('results'):
            mat = pd.read_csv(f'results/mat_{comparison}.csv', index_col=0)
            mat = mat.T
        else:
            mat = dds[[interest]].T.rename(index={interest: comparison})
            
            try:
                mat.T.to_csv(f'results/mat_{comparison}.csv')
            except OSError:
                print(f"coudn't save at results/mat_{comparison}.csv. Currently at {os.system('pwd')}")

        self.mats[comparison]= mat
        return mat
    def get_tf_df(self, mat, comparison='comparison'):
        if comparison in self.tf_dfs and self.tf_dfs[comparison] is not None:
            tf_df = self.tf_dfs[comparison]
        elif f'tf_acts_{comparison}.csv' in os.listdir('results'):
            tf_df = pd.read_csv(f'results/tf_acts_{comparison}.csv', index_col=0)
            self.tf_dfs[comparison] = tf_df
        else:
            tf_acts, tf_pvals = dc.run_ulm(mat=mat, net=self.collectri, verbose=True)
            tf_df = pd.DataFrame(tf_acts.T)
            tf_df['pvals']=tf_pvals.T
            try:
                tf_df.to_csv(f'results/tf_acts_{comparison}.csv')
            except OSError:
                print(f"coudn't save at rresults/tf_acts_{comparison}.csv. Currently at {os.system('pwd')}")
            self.tf_dfs[comparison] = tf_df
        return tf_df

    def get_tf_df_of_comparison(self, dds, interest='stat', comparison='comparison'):
        """
        This function will take a dataframe with the dds results and will return a dataframe with the DE genes
        """
        mat = self.get_mat(dds=dds, interest=interest, comparison=comparison)
        tf_df = self.get_tf_df(mat=mat, comparison=comparison)
        return tf_df
    def get_tf_df_of_combined_comparisons(self):
        all_tf_dic = {}
        for comparison, dds in self.dds_dict.items():
            tf_df = self.get_tf_df_of_comparison(dds=dds, comparison=comparison)
            tf_df = tf_df[tf_df['pvals']<self.pvalue]
            all_tf_dic[comparison] = tf_df[comparison]
            #print( tf_df.head())
        all_scores_tf_df=None
        for comparison, score in all_tf_dic.items():
            if all_scores_tf_df is None:
                all_scores_tf_df = score.rename(comparison).to_frame()
            else:
                all_scores_tf_df = all_scores_tf_df.join(score.rename(comparison))
        all_scores_tf_df.fillna(0, inplace=True)

        return all_scores_tf_df
        
        if isinstance(genes, pd.Series):  # If multiple values exist
            gene_list.append(genes.tolist())
        else:  # Single value case
            gene_list.append([genes])
    
    def save_tf_df_of_combined_comparisons(self, path='results/tf_df_combined.csv'):
        tf_df = self.get_tf_df_of_combined_comparisons()
        tf_df.to_csv(path)
        
    def get_enrtiched_tf(self, comparison='comparison', pvalue=0.1):
        if self.tf_dfs[comparison] is None:
            raise ValueError('You need to run the get_tf_df() function first')
        tf_df = self.tf_dfs[comparison]
        tf_df = tf_df[tf_df['pvals']<pvalue]
        return tf_df
    
    def get_highes_lowest_tf(self, comparison='comparison'):
        if self.tf_dfs[comparison] is None:
            raise ValueError('You need to run the get_tf_df() function first')
        tf_df = self.get_enrtiched_tf(comparison=comparison)
        tf_acts=pd.DataFrame(tf_df[comparison])
        tf_acts=tf_acts.T
        values = tf_acts.iloc[0]
        down_reg = values.sort_values(ascending=True)[:5].index.to_list()
        up_reg = values.sort_values(ascending=False)[:5].index.to_list()
        up_down_reg = down_reg.copy()
        up_down_reg.extend(up_reg)
        return up_down_reg

