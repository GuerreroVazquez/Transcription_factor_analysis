import pandas as pd
import scanpy as sc
import decoupler as dc
import matplotlib.pyplot as plt
import seaborn as sns
import textwrap

ora_fdr = 0.15
gsea_fdr=0.01
score_name='stat'
def plot_gsea_results(gsea_scores, msigdb, score_name=score_name, term=None):
    """
    Plots the top GSEA results based on NES and includes leading edge genes.
    
    Parameters:
    gsea_df: DataFrame containing GSEA results with columns Term, ES, NES, NOM p-value, FDR p-value, Set size, Tag %, Rank %, Leading edge
    msigdb_df: DataFrame with gene symbols, collections, and genesets
    top_n: Number of top terms to display
    
    Returns:
    A bar plot with the top terms and their NES.
    """
    geneset= msigdb[msigdb['geneset']==term]
    fig = dc.plot_running_score(
        df=gsea_scores,
        stat=score_name,
        net=geneset,
        source='geneset',
        target='genesymbol',
        set_name=term,
        return_fig=True
    )
    return fig

def plot_ora_results(ora_df, top_n=10, figsize=(10, 6), scale_odds_ratio=100, 
                     fontsize_title=16, fontsize_subtitle=14, fontsize_text=12,
                     fdr_bar_size=[0.8, 0.15, 0.03, 0.5], bbox_to_anchor=(1.3, 1.05), title=None):
    """
    Creates a customized lollipop plot for ORA results, with bubble size for Odds Ratio and color for FDR p-value,
    including sample circles representing minimum, maximum, and median Odds Ratios positioned above the FDR bar.
    
    Parameters:
    ora_df: DataFrame containing ORA results with columns Term, Set size, Overlap ratio, p-value, FDR p-value, Odds ratio, Combined score, Features
    top_n: Number of top terms to display
    fdr_bar_seize: array with the size and postion of the FDR bar (left, bottom, width, height)
    
    
    Returns:
    fig: The figure object for further manipulation if needed.
    ax: The axes object with the plot.
    """

    if title is None:
        title= f"Top {top_n} ORA Results by Combined Score"
    # Sort ORA results by Combined Score and select top N
    ora_top = ora_df.sort_values(by='Combined score', ascending=False).head(top_n).copy()
    
    # Replace underscores in 'Term' with spaces for better display
    ora_top['Term'] = ora_top['Term'].str.replace('_', ' ')
    
    # Wrap long labels onto two lines if necessary
    ora_top['Term'] = ora_top['Term'].apply(lambda x: '\n'.join(textwrap.wrap(x, width=30)))
    
    # Normalize FDR p-value to range between 0 and 1 for color mapping
    norm = plt.Normalize(ora_top['FDR p-value'].min(), ora_top['FDR p-value'].max())
    
    # Calculate min, max, and median Odds Ratios
    min_odds = ora_top['Odds ratio'].min()
    max_odds = ora_top['Odds ratio'].max()
    median_odds = ora_top['Odds ratio'].median()
    
    # Create the figure and axis with a fixed width
    fig, ax = plt.subplots(figsize=figsize)  # Adjust the width (10) to your preference
    
    # Bubbles for Combined Score, sized by Odds Ratio and colored by FDR p-value
    bubble = ax.scatter(ora_top['Combined score'], ora_top['Term'], 
                         s=ora_top['Odds ratio'] * scale_odds_ratio,  # Scale the size by Odds Ratio
                         c=ora_top['FDR p-value'], cmap='viridis', norm=norm, alpha=0.8)
    
    # Colorbar for FDR p-values
    cbar = plt.colorbar(bubble, ax=ax)
    cbar.set_label('FDR p-value')
    
    # Remove the black outline from the bubbles
    bubble.set_edgecolors('none')
    
    # Add labels and title, increasing the font size
    ax.set_xlabel("Combined Score", fontsize=fontsize_subtitle)
    ax.set_ylabel("Pathway", fontsize=fontsize_subtitle)
    ax.set_title(title, fontsize=fontsize_title)
    
    # Increase the font size for the y-axis labels (pathway names)
    ax.set_yticklabels(ora_top['Term'], fontsize=fontsize_text)
    
    # Add vertical gridlines
    ax.grid(True, axis='x', linestyle='--', alpha=0.7)
    
    # Remove horizontal gridlines
    ax.grid(False, axis='y')
    
    # Remove the extra square around the plot
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    
    # Add sample circles for min, median, and max Odds Ratios on top of the FDR colorbar
    sample_odds_ratios = [min_odds, max_odds]
    sample_sizes = [min_odds * scale_odds_ratio, max_odds * scale_odds_ratio]
    sample_labels = ['Min', 'Max']
    
    # Position the sample bubbles above the colorbar
    for i, (size, label, odds) in enumerate(zip(sample_sizes, sample_labels, sample_odds_ratios)):
        ax.scatter([], [], s=size, color='grey', alpha=0.5, label=f'{int(odds)}')
    
    # Legend for the sample Odds Ratio bubbles, place above the colorbar
    leg = ax.legend(scatterpoints=1, frameon=False, labelspacing=1, title='Odds ratio', 
                    loc='upper center', bbox_to_anchor=bbox_to_anchor)
    
    # Adjust the colorbar position to make room for the sample circles
    cbar.ax.set_position(fdr_bar_size)  # 
    
    # Show plot with tight layout for better fitting
    plt.tight_layout()
    
    return fig, ax