import pandas as pd

class Preprocessor:
    def __init__(self, dataframe: pd.DataFrame):
        self.dataframe = dataframe

    def calculate_m2_m1_ratios(self):
        """Calculates various m2/m1 ratios and adds them as new columns."""
        # Obtaining the data for the core ratio with the two first biopsies at the tumor core site
        self.dataframe['m2/m1 at the core'] = 0
        for i in range(1,3):
            self.dataframe['m2/m1 at the core'] += self.dataframe[f'm2_{i}'] / self.dataframe[f'm1_{i}']
        self.dataframe['m2/m1 at the core'] /= 2
        # for the Edge value, we consider the outer layers of biopsies
        self.dataframe['m2/m1 at the edge'] = 0
        for i in range(3, 6):  # Random ratios for indices 1 to 5
            self.dataframe['m2/m1 at the edge'] += self.dataframe[f'm2_{i}'] / self.dataframe[f'm1_{i}']
        self.dataframe['m2/m1 at the edge'] /= 2
        # Random aggregated ratio
        self.dataframe['random m2/m1'] = 0
        for i in range(1, 6):  # Random ratios for indices 1 to 5
            self.dataframe['random m2/m1'] += self.dataframe[f'm2_rand{i}'] / self.dataframe[f'm1_rand{i}']
        self.dataframe['random m2/m1'] /= 5

        return self.dataframe
