# This notebook contains a collection of utility functions for analyzing and visualizing logic models. 
# These functions are designed to facilitate various analysis tasks 
# including attractor analysis, model comparison, visualization, and correlation with external data.

#################################### Import libraries ####################################
# These libraries are all provided by the CoLoMoTo notebook.
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import biolqm
import ginsim
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter
from rpy2.robjects import robjects
boolnet = importr("BoolNet")
from scipy.stats import pearsonr

#################################### Attractor Analysis ####################################

def get_attractors(model_path, updating = 'asynchronous', save_path=None, model_name="model"):
    """
    Compute attractors for a Boolean network model and optionally save to CSV.
    
    Parameters:
    -----------
    model_path : str
        Path to the model file
    updating : str
        Type of updating method to use ('asynchronous' or 'synchronous')
    save_path : str, optional
        Path to save the attractor data
    model_name : str, optional
        Name of the model for file naming
        
    Returns:
    --------
    pandas.DataFrame
        DataFrame containing attractor data
    """
    # Get the attractors using BoolNet
    model = boolnet.loadNetwork(model_path)
    attr = boolnet.getAttractors(model, type=updating)
    
    # Activate the conversion context to use pandas DataFrame
    pandas2ri.activate()
    
    # Access the dataframe stored under the key '1'
    attrr = boolnet.plotAttractors(attr)
    r_df = attrr.rx2('1')
    
    # Use the local converter context to manage the conversion to a numpy array
    with localconverter(robjects.default_converter + pandas2ri.converter):
        np_array = np.array(r_df)
    
    # Extract names from the R dataframe
    row_names = list(r_df.rownames)
    column_names = list(r_df.colnames)
    
    # Create a pandas DataFrame from the numpy array
    df = pd.DataFrame(np_array, index=row_names, columns=column_names)

    if updating == 'asynchronous':
        # Add the complex/loose attractors
        attrr_comp = boolnet.plotAttractors(attr, mode='graph')
        for i in range(attrr.rx2('1').ncol, len(attrr_comp)):
            attr_array = attrr_comp[i][8][2][0]
            if len(attr_array[0]) == len(df):
                # Loop through each binary string in the array
                for idx, binary_str in enumerate(attr_array):
                    # Convert each character in the string to an integer and add as a new column
                    column_name = f'Attr{i+1}.{idx + 1}'
                    df[column_name] = [int(char) for char in binary_str]
            else:
                print("The number of digits in the binary strings does not match the number of DataFrame rows.")

    df = df.T
    df.index = df.index.str.replace('Attr', 'S')
    
    # Save to CSV if path is provided
    if save_path:
        df.to_csv(f"{save_path}/attr_{model_name}.csv")
    
    return df

#################################### Visualization ####################################

def visualize_model(model_path):
    lqm1 = biolqm.load(model_path)
    lrg1 = biolqm.to_ginsim(lqm1)
    ginsim.show(lrg1)

def plot_attractors_heatmap(df, title=None, figsize=(12, 8), cmap=None):
    """
    Create a heatmap visualization of attractors.
    
    Parameters:
    -----------
    df : pandas.DataFrame
        DataFrame with attractors
    title : str, optional
        Title for the plot
    figsize : tuple, optional
        Figure size (width, height)
    cmap : colormap, optional
        Matplotlib colormap for the heatmap
        
    Returns:
    --------
    matplotlib.figure.Figure
        The generated figure
    """
    if cmap is None:
        cmap = ListedColormap(['#c0392b', '#2980b9'])
    
    plt.figure(figsize=figsize)
    
    # Create heatmap
    ax = sns.heatmap(df, cmap=cmap, cbar=False, linewidths=0.5)
    if title:
        plt.title(title)
    
    plt.tight_layout()
    return plt.gcf()

#################################### Correlation with data ####################################

def correlation_with_data(df, data):
    """
    Calculate the Pearson correlation between the frequency of expression of a model and experimental data.
    
    Parameters:
    -----------
    df : pandas.DataFrame
        DataFrame with attractors
    data : pandas.DataFrame
        DataFrame with experimental data
        In the format of [Gene, Expression]
        
    Returns:
    --------
    plot
        A plot of the correlation
    """
    # Calculate frequency of expression for each gene
    frequency = df.sum() / len(df)
    frequency_df = frequency.reset_index()
    frequency_df.columns = ['Gene', 'Frequency']
    
    # Merge with the expression DataFrame
    merged_data = pd.merge(frequency_df, data, on='Gene', how='inner')
    
    # Compute Pearson correlation and p-value
    correlation, p_value = pearsonr(merged_data['Frequency'], merged_data['Expression'])
    
    # Plotting
    plt.figure(figsize=(6, 6))
    ax = sns.regplot(data=merged_data, x='Frequency', y='Expression', ci=None, scatter_kws={"s": 50})
    plt.xlabel('Modeled frequency of expression',fontsize = 13)
    plt.ylabel('Measured frequency of expression',fontsize = 13)
    
    # Annotate Pearson correlation and p-value
    plt.text(0.05, 0.95, f'Pearson r: {correlation:.2f}\np-value: {p_value:.2e}',
             ha='left', va='top', transform=ax.transAxes, fontsize=15, bbox=dict(boxstyle="round,pad=0.3", edgecolor='black', facecolor='white'))
    
    # Annotate gene names with adjustText to avoid overlaps
    texts = []
    for i, point in merged_data.iterrows():
        texts.append(ax.text(point['Frequency'], point['Expression'], point['Gene'], ha='right', va='center', fontsize=14))
    
    plt.show()


def calculate_phenotype_scores(df, phenotype_definitions):
    """
    Calculate phenotype scores based on activation and inhibition effects.
    
    Parameters:
    -----------
    df : pandas.DataFrame
        DataFrame with attractor states
    phenotype_definitions : dict
        Dictionary mapping phenotypes to their activators and inhibitors
        Example: {'Apoptosis': {'activators': ['TP53'], 'inhibitors': ['BCL2']}}
        
    Returns:
    --------
    pandas.DataFrame
        DataFrame with phenotype scores
    """
    scores = pd.DataFrame(index=df.index)
    
    for phenotype, definition in phenotype_definitions.items():
        activators = definition.get('activators', [])
        inhibitors = definition.get('inhibitors', [])
        
        # Calculate score for each state/attractor
        score = pd.Series(0, index=df.index)
        
        for activator in activators:
            if activator in df.columns:
                score += df[activator]
                
        for inhibitor in inhibitors:
            if inhibitor in df.columns:
                score -= df[inhibitor]
                
        scores[phenotype] = score
    
    return scores

def correlate_with_outcome(phenotype_scores, outcome_data, outcome_col='Outcome'):
    """
    Correlate phenotype scores with clinical outcomes
    
    Parameters:
    -----------
    phenotype_scores : pandas.DataFrame
        DataFrame with phenotype scores for each attractor
    clinical_data : pandas.DataFrame
        DataFrame with clinical data
    outcome_col : str, optional
        Column name for outcome variable
        
    Returns:
    --------
    pandas.DataFrame
        DataFrame with correlation results
    """
    
    # Initialize results dataframe
    results = []
    
    # Loop through phenotypes
    for phenotype in phenotype_scores.columns:
        # Compute correlation between phenotype score and outcome
        corr, p_value = pearsonr(
            phenotype_scores[phenotype],
            outcome_data[outcome_col]
        )
        
        results.append({
            'Phenotype': phenotype,
            'Correlation': corr,
            'P_value': p_value
        })
    
    return pd.DataFrame(results)

#################################### model conversion ####################################
def convert_txt_to_sbml(txt_file_path, sbml_file_path):
    """
    Convert Boolean network in text format to SBML format.
    
    Parameters:
    -----------
    txt_file_path : str
        Path to the text file
    sbml_file_path : str
        Path for saving the SBML file
    
    Returns:
    --------
    None
    """
    # Load the network from text file
    network = boolnet.loadNetwork(txt_file_path)
    
    # Convert to SBML format
    boolnet.toSBML(network, sbml_file_path)
    
    print(f"Converted {txt_file_path} to {sbml_file_path}")
    return

def convert_sbml_to_txt(sbml_file_path, txt_file_path):
    """
    Convert Boolean network in SBML format to text format.
    """
    # Load the network from SBML file
    network = boolnet.loadSBML(sbml_file_path)
    boolnet.saveNetwork(network, txt_file_path)
    print(f"Converted {sbml_file_path} to {txt_file_path}")
    return