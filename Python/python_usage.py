# Example Usage
data = pd.DataFrame({
    "Date": pd.date_range(start="2020-01-01", periods=208, freq="W"),
    "EQ_Volume": np.random.normal(100, 20, 208),
    "Media_1": np.random.normal(50, 10, 208),
    "Media_2": np.random.normal(30, 5, 208)
})

dates = data["Date"]
distribution_params = {
    "Weibull": {"param1": [0.1, 5, 0.1], "param2": [1, 70, 1]},
    "Lognormal": {"param1": [1, 20, 1], "param2": [10, 20, 1]},
    "Logistic": {"param1": [1, 16, 1], "param2": [1, 16, 1]},
    "NegExp": {"param1": [0.1, 1, 0.1], "param2": [None, None, None]},
    "ChiSquare": {"param1": [1, 70, 1], "param2": [None, None, None]},
}

saturation_models = {
    "linear": lambda x: x,
    "diminishing": lambda x: x / (1 + x),
    "log": lambda x: np.log1p(x),
    "exponential": lambda x: 1 - np.exp(-x)
}

adstock_model = AdstockModel(data, dates, saturation_models, distribution_params)
adstock_model.run_adstock()
adstock_model.results.to_csv("adstock_new_iterations.csv", index=False)
