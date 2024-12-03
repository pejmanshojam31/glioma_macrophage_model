
import pandas as pd

class DataLoader:
    def __init__(self, file_path: str):
        self.file_path = file_path

    def load_data(self):
        try:
            df = pd.read_csv(self.file_path)
            print(f"Dataset loaded successfully from {self.file_path}")
            return df
        except Exception as e:
            print(f"Error loading dataset: {e}")
            raise
        