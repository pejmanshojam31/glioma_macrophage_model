'''
This code is used to obtain Figures 6 and 7 of the manuscript.
 We first need to analyse the correlation coefficient values in a heatmap
between the features and the target variable.
 We then calculate the mutual information to understand more about the nature of our dataset.
'''


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.feature_selection import mutual_info_regression

class DataAnalyzer:
    def __init__(self,dataframe: pd.DataFrame):
        self.dataframe = dataframe
    def calculate_correlation_matrix(self,save_path=None):
        """
        Calculates and visualizes the correlation matrix for the dataset.

        Parameters:
            save_path (str, optional): Path to save the heatmap as a file. Defaults to None.

        Returns:
            pd.DataFrame: Correlation matrix.
        """
        correlation_matrix=self.dataframe.corr()

        # Plotting the heatmap
        plt.figure(figsize=(12,10))
        sns.heatmap(correlation_matrix, annot=True, fmt=".2f", cmap='Blues', 
                    annot_kws={'size':12, 'weight': 'bold'})
        plt.title('Correlation Matrix', fontsize=18, fontweight='bold')
        plt.xticks(fontsize=14,fontweight='bold')
        plt.yticks(fontsize=14,fontweight='bold')

        # Save to file if save_path is provided
        if save_path:
            plt.savefig(save_path,format='pdf')
        plt.show()

        return correlation_matrix

    def calculate_and_plot_mutual_info(self, features, target, target_name, color='teal', save_path=None):
        """
        Calculates and visualizes mutual information between features and a target variable.

        Parameters:
            features (DataFrame): Feature set.
            target (array-like): Target variable.
            target_name (str): Name of the target variable for plot title.
            color (str): Color for the bar plot.
            save_path (str, optional): File path to save the figure. Defaults to None.
        """
        # Calculate mutual information
        mutual_info = mutual_info_regression(features, target)
        mutual_info_series = pd.Series(mutual_info, index=features.columns)
        sorted_mutual_info = mutual_info_series.sort_values(ascending=False)

        # Plot mutual information
        plt.figure(figsize=(10, 6))
        sorted_mutual_info.plot(kind='bar', color=color)
        plt.title(f'Mutual Information with {target_name}', fontsize=16, fontweight='bold')
        plt.ylabel('Mutual Information Value', fontsize=14, fontweight='bold')
        plt.xticks(rotation=45, fontsize=12, fontweight='bold')
        plt.yticks(fontsize=12, fontweight='bold')
        plt.grid(axis='y', linestyle='--', alpha=0.7)
        plt.tight_layout()

        # Save the figure if save_path is provided
        if save_path:
            plt.savefig(save_path, format='pdf')
            print(f"Figure saved to {save_path}")

        # Show the plot
        plt.show()