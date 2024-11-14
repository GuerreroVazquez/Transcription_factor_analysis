import pandas as pd
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib import rcParams
import seaborn as sns
from scipy.stats import ttest_ind



    # Define age groups
bins = [18, 35, 65, 100]
labels = ['Young', 'Middle Age', 'Old']
palette = {'Male': 'blue', 'Female': 'pink'}
sig_p_value = 0.05
default_color = 'green'
def get_gene_and_age(df, gene):
    gene_data = df[['Age', gene]].dropna()
    return gene_data

def plot_expression_over_age(gene_data, gene,save=None, plot=True):
    plt.figure(figsize=(8, 6))
    plt.plot(gene_data['Age'], gene_data[gene], marker='o')
    plt.title(f"Expression of {gene} over Age")
    plt.xlabel("Age")
    plt.ylabel(f"Expression of {gene}")
    plt.grid(True)
    if save:
        plt.savefig(save)
    plt.show()

def scatter_plot_expression_over_age(gene_data, gene, save=None, plot=True):
    plt.figure(figsize=(8, 6))
    sns.scatterplot(x=gene_data['Age'], y=gene_data[gene], hue=gene_data['Age'], palette='viridis', legend=False)
    plt.title(f"Expression of {gene} across Age")
    plt.xlabel("Age")
    plt.ylabel(f"Expression of {gene}")
    plt.grid(True)
    if save:
        plt.savefig(save)
    plt.show()

def violin_plot_grouped_by_age(gene_data, gene, save=None, plot=True, color='green'):
    
    gene_data.loc[:, 'Age Group'] = pd.cut(gene_data['Age'], bins=bins, labels=labels, right=False)
    
    # Plot violin plot
    fig, ax = plt.subplots(figsize=(8, 6))
    #plt.figure(figsize=(8, 6))
    sns.violinplot(x='Age Group', y=gene, data=gene_data, color=color)
    ax.set_title(f"Violin plot of {gene} grouped by Age")
    ax.set_xlabel("Age Group")
    ax.set_ylabel(f"Expression of {gene}")
    ax.grid(True)
    if save:
        plt.savefig(save)
    if plot:
        plt.show()
    return fig

def gene_profile(df, gene, save=None, plot=True):
    save_sp = None
    save_vp = None
    if save:
        save_sp = save + "_sp.jpg"
        save_vp = save + "_vp.jpg"
        
    gene_data = get_gene_and_age(df, gene)  # Extract gene data
    scatter_plot_expression_over_age(gene_data, gene, save_sp)  # Plot expression over age
    violin_plot_grouped_by_age(gene_data, gene, save_vp)  # Make violin plots


# Function 4: Violin plot grouped by sex
def violin_plot_grouped_by_sex(gene_data, gene, save = None, plot=True, color=None):

    fig, ax = plt.subplots(figsize=(8, 6))
    # Plot violin plot
    #plt.figure(figsize=(8, 6))
    sns.violinplot(x='Sex', y=gene, data=gene_data, palette=palette)
    ax.set_title(f"Violin plot of {gene} grouped by Sex")
    ax.set_xlabel("Sex")
    ax.set_ylabel(f"Expression of {gene}")
    ax.grid(True)
    if save:
        plt.savefig(save)
    if plot:
        plt.show()
    return fig

# Function 5: Violin plot grouped by experiment
def violin_plot_grouped_by_experiment(gene_data, gene, save = None, plot=True, color='green'):
    # Plot violin plot
    fig, ax = plt.subplots(figsize=(8, 6))
    #plt.figure(figsize=(8, 6))
    sns.violinplot(x='Experiment', y=gene, data=gene_data, color=color)
    ax.set_title(f"Violin plot of {gene} grouped by Experiment")
    ax.set_xlabel("Experiment")
    ax.set_ylabel(f"Expression of {gene}")
    ax.grid(True)
    if save:
        plt.savefig(save)
    if plot:
        plt.show()
    return fig

# Function 6: Violin plot grouped by sex and age group
def violin_plot_grouped_by_sex_and_age_group(gene_data, gene, save=None, plot=True, color=None):
    


    gene_data['Age Group'] = pd.cut(gene_data['Age'], bins=bins, labels=labels, right=False)
    
    # Plot violin plot
    fig, ax = plt.subplots(figsize=(8, 6))
    #plt.figure(figsize=(8, 6))
    sns.violinplot(x='Age Group', y=gene, hue='Sex', data=gene_data, split=True, palette=palette)
    ax.set_title(f"Violin plot of {gene} grouped by Sex and Age Group")
    ax.set_xlabel("Age Group")
    ax.set_ylabel(f"Expression of {gene}")
    ax.grid(True)
    if save:
        plt.savefig(save)
    if plot:
        plt.show()
    return fig



def calculate_pvalue_gene_change(gene_data, gene):
    # Use .loc to modify the DataFrame safely
    gene_data.loc[:, 'Age Group'] = pd.cut(gene_data['Age'], bins=bins, labels=labels, right=False)
    
    # Split data into age groups
    young = gene_data[gene_data['Age Group'] == 'Young'][gene]
    middle_age = gene_data[gene_data['Age Group'] == 'Middle Age'][gene]
    old = gene_data[gene_data['Age Group'] == 'Old'][gene]
    
    # Dictionary to store p-values
    pvalues = {}
    
    # Perform T-tests
    pvalues['Young_vs_MiddleAge'] = ttest_ind(young, middle_age, nan_policy='omit').pvalue
    pvalues['MiddleAge_vs_Old'] = ttest_ind(middle_age, old, nan_policy='omit').pvalue
    pvalues['Young_vs_Old'] = ttest_ind(young, old, nan_policy='omit').pvalue
    
    return pvalues




def generate_three_letter_code(gene_data, gene, pvalues):
    # Define age groups again (to calculate means)

    gene_data.loc[:, 'Age Group'] = pd.cut(gene_data['Age'], bins=bins, labels=labels, right=False)
    
    # Split data into age groups
    young = gene_data[gene_data['Age Group'] == 'Young'][gene]
    middle_age = gene_data[gene_data['Age Group'] == 'Middle Age'][gene]
    old = gene_data[gene_data['Age Group'] == 'Old'][gene]
    
    def compare_means(group1, group2):
        # Compare means of two groups
        if group1.mean() < group2.mean():
            return 'B'  # Younger group's mean is smaller
        else:
            return 'C'  # Younger group's mean is larger or equal
    
    # Initialize the 3-letter code
    code = ""
    
    # Compare Young vs. Middle Age
    if pvalues['Young_vs_MiddleAge'] >= sig_p_value:
        code += 'A'  # Not significant
    else:
        code += compare_means(young, middle_age)
    
    # Compare Middle Age vs. Old
    if pvalues['MiddleAge_vs_Old'] >= sig_p_value:
        code += 'A'
    else:
        code += compare_means(middle_age, old)
    
    # Compare Young vs. Old
    if pvalues['Young_vs_Old'] >= sig_p_value:
        code += 'A'
    else:
        code += compare_means(young, old)
    
    return code


def get_change_curve(gene_data, gene):
    pvals = calculate_pvalue_gene_change(gene_data, gene)
    code = generate_three_letter_code(gene_data, gene, pvals)
    return code


def load_gene_data_with_metadata(csv_file):
    # Read the CSV file without header (since the first 3 rows are metadata)
    df = pd.read_csv(csv_file, header=None)
    
    # Dictionary to store the data
    gene_data = {}

    # Iterate over each column
    for col in df.columns:
        # Extract metadata from the first 3 rows
        metadata = {
            'Algorithm': df.iloc[0, col],  # First row as metadata 1
            'Extra': df.iloc[1, col],  # Second row as metadata 2
            'Selection': df.iloc[2, col],  # Third row as metadata 3
        }
        
        # Extract the gene data from the 4th row onwards
        genes = df.iloc[3:, col].dropna().tolist()  # Drop any NaNs (for varying lengths)
        
        # Store the metadata and genes in the dictionary
        gene_data[col] = {
            'metadata': metadata,
            'genes': genes
        }

    return gene_data

def create_violin_pdf(gene_data, genes, output_file, title="Gene Violin Plots", violin_function=violin_plot_grouped_by_age, color='green'):
    with PdfPages(output_file) as pdf:
        # Create the index page
        fig, ax = plt.subplots(figsize=(8.5, 11))  # Letter size page
        
        ax.set_frame_on(False)
        ax.set_xticks([])
        ax.set_yticks([])
        
        # Title for the index
        ax.text(0.5, 0.95, title, ha='center', fontsize=16, fontweight='bold')
        
        # Generate index with clickable links (not yet, just a visual index)
        y_position = 0.9
        for i, gene in enumerate(genes):
            # Add gene name to the index
            ax.text(0.1, y_position - (i * 0.04), f"{i+1}. {gene}", fontsize=12, ha='left')
        
        # Save index page
        pdf.savefig(fig)
        plt.close(fig)

        # Create violin plots for each gene using the modified function
        for gene in genes:
            # Generate the violin plot and return the figure
            fig = violin_function(gene_data, gene, color=color)
            
            # Save the figure to the PDF
            pdf.savefig(fig)
            
            # Close the figure after saving it to free up memory
            plt.close(fig)
