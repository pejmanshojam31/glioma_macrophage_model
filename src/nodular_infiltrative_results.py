import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import mean_squared_error
from joblib import Parallel, delayed
from sklearn.ensemble import GradientBoostingRegressor

class NodularAndInfiltrativeResults:
    def __init__(self, dataframe, output_path, random_seed=42):
        """
        Initialize the class with the main DataFrame and output path.

        Parameters:
            dataframe (pd.DataFrame): Preprocessed dataset.
            output_path (str): Path to save results.
        """
        self.dataframe = dataframe
        self.output_path = output_path
        self.feature_importances = {}
        self.random_seed = random_seed
    def split_datasets(self):
        """Split the dataset into nodular and infiltrative subsets based on the 25th percentile of 'IW'."""
        Q1 = self.dataframe['IW'].quantile(0.25)
        self.nodular = self.dataframe[self.dataframe['IW'] <= Q1].copy()
        self.infiltrative = self.dataframe[self.dataframe['IW'] > Q1].copy()

    def add_noise_to_datasets(self, features, noise_levels, verbose=False):
        """
        Generate datasets with varying noise levels for both nodular and infiltrative subsets.

        Parameters:
            features (list): Features to which noise will be added.
            noise_levels (list): List of noise levels (sigma values).
            verbose (bool): Whether to print processing details. Default is False.
        """
        def add_noise(dataset, sigma):
            noisy_dataset = dataset.copy()
            rng = np.random.default_rng(seed=self.random_seed)  # Use a reproducible random generator
            for feature in features:
                if feature in dataset.columns:
                    noise = rng.normal(0, sigma, size=dataset[feature].shape)
                    noisy_dataset[feature] += noise
                else:
                    if verbose:
                        print(f"Feature {feature} not found in dataset. Skipping.")
            return noisy_dataset

        # Initialize the datasets with "No Noise" (original copies)
        self.datasets_nodular = {"No Noise": self.nodular.copy()}
        self.datasets_infiltrative = {"No Noise": self.infiltrative.copy()}

        # Generate noisy datasets for each specified noise level
        for sigma in noise_levels:
            if sigma > 0:  # Skip noise generation for "No Noise"
                if verbose:
                    print(f"Adding noise with sigma={sigma} to nodular and infiltrative datasets.")
                self.datasets_nodular[f"{sigma} Noise"] = add_noise(self.nodular, sigma)
                self.datasets_infiltrative[f"{sigma} Noise"] = add_noise(self.infiltrative, sigma)

    def run_analysis(self,scenarios,targets):
        """
        Perform the analysis for nodular and infiltrative datasets.

        Parameters:
            scenarios (dict): Dictionary of feature scenarios.
            targets (list): List of target columns (e.g., ['IW_at1', 'IW_at3']).

        Returns:
            pd.DataFrame: Results of the analysis.
        """
        self.scenarios=scenarios
        self.all_datasets ={
            "Nodular":self.datasets_nodular,
            "Infiltrative":self.datasets_infiltrative
        }

        results = [
            result for result in Parallel(n_jobs=-1)(
                delayed(self._process_dataset)(
                    dataset_type, dataset_name,dataset,scenario_name, features, target
                )
                for dataset_type, datasets in self.all_datasets.items()
                for dataset_name, dataset in datasets.items()
                for scenario_name, features in scenarios.items()
                for target in targets
            )
            if result  # Exclude empty results
        ]

        if not results:
            print("No valid results generated. Check your scenarios, targets, and datasets.")
            return pd.DataFrame()

        results_df = pd.DataFrame(results)
        results_df.to_csv(self.output_path,index=False)
        print(f"Results saved to {self.output_path}")
        return results_df


    def _process_dataset(self, dataset_type, dataset_name, dataset, scenario_name, features, target):
        """
        Helper method to process a single dataset for analysis.
        Parameters:
            dataset_type (str): Type of dataset ('Nodular' or 'Infiltrative').
            dataset_name (str): Name of the dataset.
            dataset (pd.DataFrame): Dataset to process.
            scenario_name (str): Name of the scenario.
            features (list): Features for the scenario.
            target (str): Target variable.

        Returns:
            dict: A dictionary containing the results.
        """
        # Initialize the result with default values
        result = {
            "Type": dataset_type,
            "Scenario": scenario_name,
            "Dataset": dataset_name,
            "TimePoint": target,
            "MSE": None,
            "Status": "Failed"
        }

        try:
            # Check for missing features or target
            missing_features = [feature for feature in features if feature not in dataset.columns]
            if missing_features:
                print(f"Skipping scenario '{scenario_name}' for dataset '{dataset_name}' due to missing features: {missing_features}")
                result["Status"] = f"Missing features: {missing_features}"
                return result
            
            if target not in dataset.columns:
                print(f"Skipping target '{target}' for dataset '{dataset_name}' and scenario '{scenario_name}' due to missing target column.")
                result["Status"] = f"Missing target: {target}"
                return result

            # Extract features and target
            X = dataset[features].values
            Y = dataset[target].values

            if X.size == 0 or Y.size == 0:
                print(f"Empty data encountered for scenario '{scenario_name}' in dataset '{dataset_name}'.")
                result["Status"] = "Empty data"
                return result

            # Split data into training and testing sets
            X_train, X_test, Y_train, Y_test = train_test_split(
                X, Y, test_size=0.2, random_state=self.random_seed
            )

            # Standardize features
            scaler = StandardScaler()
            X_train_scaled = scaler.fit_transform(X_train)
            X_test_scaled = scaler.transform(X_test)

            # Train model and calculate MSE
            gb_optimized = GradientBoostingRegressor(
                learning_rate=0.05, max_depth=3, n_estimators=200, random_state=self.random_seed
            )
            gb_optimized.fit(X_train_scaled, Y_train)
            Y_pred = gb_optimized.predict(X_test_scaled)
            mse = mean_squared_error(Y_test, Y_pred)

            # Populate the result dictionary with successful results
            result["MSE"] = mse
            result["Status"] = "Success"

        except Exception as e:
            print(f"An error occurred while processing scenario '{scenario_name}' for dataset '{dataset_name}': {e}")
            result["Status"] = f"Error: {e}"

        return result


    def plot_mse(self, results_df, scenarios, dataset_type, output_prefix):
        """
        Plot MSE for the specified dataset type using the Viridis colormap.

        Parameters:
            results_df (pd.DataFrame): DataFrame containing the results.
            scenarios (dict): Dictionary of scenarios.
            dataset_type (str): Dataset type ('Nodular' or 'Infiltrative').
            output_prefix (str): Prefix for the output file.
        """
        try:
            # Filter results for the given dataset type
            dataset_results = results_df[results_df['Type'] == dataset_type]

            # Define bar positions and unique datasets
            ind = np.arange(len(scenarios))  # Scenario positions
            width = 0.2  # Bar width
            datasets = dataset_results['Dataset'].unique()  # Unique datasets (e.g., noise levels)

            # Use Viridis colormap for bars
            mse_colors = plt.cm.viridis(np.linspace(0, 1, len(datasets)))

            # Find the minimum and maximum MSE for this dataset type
            mse_min = dataset_results['MSE'].min()
            mse_max = dataset_results['MSE'].max()

            for time_point in dataset_results['TimePoint'].unique():
                fig, ax = plt.subplots(figsize=(14, 8))

                for i, (dataset, color) in enumerate(zip(datasets, mse_colors)):
                    # Filter data for the current dataset, time point, and dataset type
                    mse_values = dataset_results[
                        (dataset_results['Dataset'] == dataset) &
                        (dataset_results['TimePoint'] == time_point)
                    ]['MSE']

                    # Ensure there are valid MSE values to plot
                    if mse_values.empty:
                        print(f"Warning: No MSE values found for {dataset} at {time_point}.")
                        continue

                    # Plot bars
                    bar_positions = ind + i * width
                    bars = ax.bar(bar_positions, mse_values, width, label=f"{dataset}", color=color)

                    # Annotate bars with MSE values
                    for bar, value in zip(bars, mse_values):
                        ax.text(
                            bar.get_x() + bar.get_width() / 2,
                            bar.get_height(),
                            f'{value:.2f}',
                            ha='center',
                            va='bottom',
                            fontsize=10,
                            fontweight='bold'
                        )

                # Configure plot
                ax.set_ylabel('MSE', fontsize=20, fontweight='bold')
                ax.set_xlabel('Scenarios', fontsize=22, fontweight='bold')
                ax.set_xticks(ind + width * (len(datasets) - 1) / 2)
                ax.set_xticklabels(scenarios.keys(), rotation=45, fontsize=20, fontweight='bold')
                
                # Set dynamic y-axis limits for the dataset type
                ax.set_ylim([mse_min, mse_max + 2])
                ax.tick_params(axis='y', labelsize=18)
                
                # Save and show plot
                fig.savefig(f"{output_prefix}_{time_point}_{dataset_type}.png", format='png', dpi=300, bbox_inches='tight')
                plt.tight_layout()
                plt.show()

        except Exception as e:
            print(f"An error occurred while plotting MSE: {e}")



    def _get_features(self, dataset):
        """
        Extract the features to be used for analysis by excluding specific columns.

        Parameters:
            dataset (pd.DataFrame): Dataset from which to extract features.

        Returns:
            list: List of feature column names.
        """
        exclude_columns = [
            "param_b", "param_d", 'param_h2', 'param_r', 'param_h3', 'param_Dm', 'param_ts', 'param_S', 'IW',
            "IW_after_3_months", "IW_after_6_months", "IW_after_12_months", 'm1_1', 'm1_2', 'm1_3', 'm1_4', 'm1_5',
            'm2_1', 'm2_2', 'm2_3', 'm2_4', 'm2_5', 'm1_rand1', 'm1_rand2', 'm1_rand3', 'm1_rand4', 'm1_rand5',
            'm2_rand1', 'm2_rand2', 'm2_rand3', 'm2_rand4', 'm2_rand5',
        ]
        return dataset.drop(columns=exclude_columns).columns.tolist()

    def calculate_feature_importances(self, targets):
        """
        Calculate feature importances for all datasets for the specified targets.

        Parameters:
            targets (list): List of target columns (e.g., ['IW_after_3_months', 'IW_after_12_months']).
        """
        self.all_datasets = {
            "Nodular": self.datasets_nodular,
            "Infiltrative": self.datasets_infiltrative
        }

        if not self.all_datasets:
            raise ValueError("No datasets available. Ensure datasets have been created using `add_noise_to_datasets`.")

        for target in targets:
            for dataset_type, datasets in self.all_datasets.items():
                for dataset_name, dataset in datasets.items():
                    features = self._get_features(dataset)

                    if target not in dataset.columns:
                        print(f"Skipping target {target} for {dataset_name} due to missing target column.")
                        continue

                    X = dataset[features].values
                    Y = dataset[target].values

                    scaler = StandardScaler()
                    X_scaled = scaler.fit_transform(X)

                    gb_optimized = GradientBoostingRegressor(learning_rate=0.05, max_depth=3, n_estimators=200, random_state=0)
                    gb_optimized.fit(X_scaled, Y)

                    if target not in self.feature_importances:
                        self.feature_importances[target] = {}
                    if dataset_type not in self.feature_importances[target]:
                        self.feature_importances[target][dataset_type] = {}
                    self.feature_importances[target][dataset_type][dataset_name] = gb_optimized.feature_importances_

        print("Feature importance calculation complete.")

    def plot_feature_importances(self, target, dataset_type, output_prefix):
        """
        Plot feature importances for the specified target and dataset type.

        Parameters:
            target (str): Target variable (e.g., 'IW_after_3_months').
            dataset_type (str): Dataset type ('Nodular' or 'Infiltrative').
            output_prefix (str): Prefix for saving the output file.
        """
        if target not in self.feature_importances:
            print(f"No feature importance data available for target: {target}")
            return
        if dataset_type not in self.feature_importances[target]:
            print(f"No feature importance data available for dataset type: {dataset_type}")
            return

        all_importances = self.feature_importances[target][dataset_type]

        sample_dataset = next(iter(self.all_datasets[dataset_type].values()))
        features = self._get_features(sample_dataset)

        positions = np.arange(len(features))
        bar_width = 0.2

        plt.figure(figsize=(12, 8))
        for idx, (noise_level, importances) in enumerate(all_importances.items()):
            if len(importances) != len(features):
                print(f"Mismatch between importances and features for {noise_level}. Skipping this noise level.")
                continue

            bar_positions = positions + idx * bar_width
            plt.bar(bar_positions, importances, bar_width, label=noise_level)

            for x, y in zip(bar_positions, importances):
                plt.text(x, y + 0.01, f'{y:.2f}', ha='center', va='bottom', fontsize=10, fontweight='bold')

        plt.ylabel('Feature Importance', fontsize=20, fontweight='bold')
        plt.xticks([p + bar_width / 2 * (len(all_importances) - 1) for p in positions],
                   features, rotation='vertical', fontsize=18, fontweight='bold')
        plt.yticks(fontsize=20, fontweight='bold')
        plt.title(f'Feature Importance for {dataset_type} ({target})', fontsize=18, fontweight='bold')
        plt.legend(title='Noise Levels', fontsize=15, title_fontsize=16, loc='upper right')
        plt.tight_layout()
        plt.savefig(f"{output_prefix}_{target}_{dataset_type}.png", format='png', dpi=300)
        plt.show()


