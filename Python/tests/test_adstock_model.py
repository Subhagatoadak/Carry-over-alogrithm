import pytest
import numpy as np
import pandas as pd
from adstock_model import AdstockModel

@pytest.fixture
def sample_data():
    data = pd.DataFrame({
        "Date": pd.date_range(start="2020-01-01", periods=10, freq="W"),
        "EQ_Volume": np.random.normal(100, 20, 10),
        "Media_1": np.random.normal(50, 10, 10),
        "Media_2": np.random.normal(30, 5, 10)
    })
    return data

@pytest.fixture
def sample_saturation_models():
    return {
        "linear": lambda x: x,
        "diminishing": lambda x: x / (1 + x),
        "log": lambda x: np.log1p(x),
        "exponential": lambda x: 1 - np.exp(-x)
    }

@pytest.fixture
def sample_distribution_params():
    return {
        "Weibull": {"param1": [0.1, 5, 0.1], "param2": [1, 70, 1]},
        "Lognormal": {"param1": [1, 20, 1], "param2": [10, 20, 1]},
        "Logistic": {"param1": [1, 16, 1], "param2": [1, 16, 1]},
        "NegExp": {"param1": [0.1, 1, 0.1], "param2": [None, None, None]},
        "ChiSquare": {"param1": [1, 70, 1], "param2": [None, None, None]}
    }

@pytest.fixture
def adstock_model(sample_data, sample_saturation_models, sample_distribution_params):
    dates = sample_data["Date"]
    return AdstockModel(sample_data, dates, sample_saturation_models, sample_distribution_params)

def test_adstock_initialization(adstock_model):
    assert isinstance(adstock_model, AdstockModel)
    assert adstock_model.data is not None
    assert isinstance(adstock_model.saturation_models, dict)
    assert isinstance(adstock_model.distribution_params, dict)

def test_adstock_run(adstock_model):
    adstock_model.run_adstock()
    assert adstock_model.results is not None
    assert not adstock_model.results.empty

def test_compute_correlations(adstock_model):
    adstock_model.run_adstock()
    correlations = adstock_model.compute_correlations()
    assert isinstance(correlations, pd.DataFrame)
    assert "Variable" in correlations.columns
    assert "Correl" in correlations.columns
    assert "Pval" in correlations.columns

def test_lagpad_function(adstock_model):
    data = np.array([1, 2, 3, 4, 5])
    lagged_data = adstock_model.lagpad(data, 2)
    assert np.isnan(lagged_data[0]) and np.isnan(lagged_data[1])
    assert lagged_data[2] == 1

def test_normalization(adstock_model):
    data = np.array([1, 2, 3, 4, 5])
    normalized_data = adstock_model.normalize(data)
    assert np.isclose(np.sum(normalized_data), 1)

def test_ztran(adstock_model):
    data = np.array([1, 2, 3, 4, 5])
    transformed_data = adstock_model.ztran(data)
    assert np.isclose(np.mean(transformed_data), 0)
    assert np.isclose(np.std(transformed_data), 1)

if __name__ == "__main__":
    pytest.main()
