'''
Here we test multiple models to find the best one for each dataset. This is preliminary and we used it to understand the model
performances.
'''
from sklearn.svm import SVR
from sklearn.neighbors import KNeighborsRegressor
from sklearn.linear_model import Ridge, ElasticNet, LinearRegression, Lasso
from sklearn.tree import DecisionTreeRegressor
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from xgboost import XGBRegressor
from sklearn.model_selection import train_test_split, cross_val_score
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import Pipeline
from sklearn.metrics import mean_squared_error
import numpy as np

class ModelSelector:
    def __init__(self, X, Y):
        """
        Initialize ModelSelector with feature and target datasets
        Parameters:
            X (DataFrame): Features dataset.
            Y (array-like): Target variable.
        """
        self.X = X
        self.Y = Y
        self.models = {
            "Linear Regression":LinearRegression(),
            "Lasso":Lasso(random_state=0),
            "Decision Tree":DecisionTreeRegressor(random_state=0),
            "Random Forest":RandomForestRegressor(random_state=0),
            "Gradient Boosting":GradientBoostingRegressor(random_state=0),
            "SVR":SVR(),
            "KNN": KNeighborsRegressor(),
            "Elastic Net":ElasticNet(random_state=0),
            "XGBoost": XGBRegressor(random_state=0),
        }
        self.cv_scores ={}
        self.mse_train_scores ={}
        self.mse_test_scores ={}

    def perform_cross_validation(self, cv=10):
        """
        Perform cross-validation for all models and store the mean MSE scores.
        
        Parameters:
            cv (int): Number of folds for cross-validation. Defaults to 10.
        """
        for name, model in self.models.items():
            pipeline = Pipeline([
                ('scaler', StandardScaler()),  # Standardization
                ('model', model)  # Model
            ])
            scores = -cross_val_score(pipeline, self.X, self.Y, scoring='neg_mean_squared_error', cv=cv)
            self.cv_scores[name] = scores.mean()

    def evaluate_models(self, test_size=0.2, random_state=0):
        """
        Train and evaluate all models, calculating MSE for both training and test sets.

        Parameters:
            test_size (float): Proportion of the dataset to include in the test split. Defaults to 0.2.
            random_state (int): Seed for reproducibility. Defaults to 0.
        """
        # Split the data
        X_train, X_test, Y_train, Y_test = train_test_split(self.X, self.Y, test_size=test_size, random_state=random_state)
        
        # Standardize the data
        scaler = StandardScaler()
        X_train_standardized = scaler.fit_transform(X_train)
        X_test_standardized = scaler.transform(X_test)

        for name, model in self.models.items():
            model.fit(X_train_standardized, Y_train)
            Y_train_pred = model.predict(X_train_standardized)
            Y_test_pred = model.predict(X_test_standardized)
            
            mse_train = mean_squared_error(Y_train, Y_train_pred)
            mse_test = mean_squared_error(Y_test, Y_test_pred)
            
            self.mse_train_scores[name] = mse_train
            self.mse_test_scores[name] = mse_test

    def display_results(self):
        """
        Display cross-validation MSE and train/test MSE for each model.
        """
        print("Cross-Validation MSE Scores:")
        for model_name, mse in self.cv_scores.items():
            print(f"{model_name}: {mse:.4f}")
        
        print("\nTraining MSE Scores:")
        for model_name, mse in self.mse_train_scores.items():
            print(f"{model_name}: {mse:.4f}")
        
        print("\nTest MSE Scores:")
        for model_name, mse in self.mse_test_scores.items():
            print(f"{model_name}: {mse:.4f}")
