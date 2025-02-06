import numpy as np
import pandas as pd
from scipy.stats import weibull_min, lognorm, logistic, expon, chi2
from scipy.stats import pearsonr

class AdstockModel:
    def __init__(self, data, dates, saturation_models, distribution_params):
        self.data = data
        self.results = None
        self.dates = dates
        self.saturation_models = saturation_models
        self.distribution_params = distribution_params
    
    def lagpad(self, x, k):
        if k > 0:
            return np.concatenate((np.full(k, np.nan), x))[:len(x)]
        else:
            return np.concatenate((x[-k:], np.full(-k, np.nan)))
    
    def normalize(self, x):
        return x / np.nansum(x) if np.nansum(x) != 0 else x
    
    def ztran(self, x):
        return (x - np.nanmean(x)) / np.nanstd(x) if np.nanstd(x) != 0 else x
    
    def generate_distributions(self, dist_name, params, func):
        t = pd.DataFrame({"X": np.arange(1, 105)})
        for i in np.arange(*params['param1']):
            for j in np.arange(*params['param2']):
                t[f"{dist_name}_{i}_{j}"] = func(t["X"], i, j)
        return t
    
    def compute_correlations(self):
        df = []
        for col in self.data.columns[1:]:
            corr, pval = pearsonr(self.data[col].dropna(), self.data.iloc[:, 0].dropna())
            df.append({"Variable": col, "Correl": corr, "Pval": pval})
        df = pd.DataFrame(df).dropna()
        return df.sort_values(by="Correl", ascending=False)
    
    def run_adstock(self):
        distributions = {
            "Weibull": self.generate_distributions("Weibull", self.distribution_params["Weibull"], weibull_min.pdf),
            "Lognormal": self.generate_distributions("Lognormal", self.distribution_params["Lognormal"], lognorm.pdf),
            "Logistic": self.generate_distributions("Logistic", self.distribution_params["Logistic"], logistic.pdf),
            "NegExp": self.generate_distributions("NegExp", self.distribution_params["NegExp"], expon.pdf),
            "ChiSquare": self.generate_distributions("ChiSquare", self.distribution_params["ChiSquare"], chi2.pdf),
        }
        
        final_results = []
        for dist_name, dist_data in distributions.items():
            dist_data.iloc[:, 1:] = dist_data.iloc[:, 1:].apply(self.normalize, axis=0)
            dist_data = pd.concat([dist_data, pd.DataFrame(np.zeros((104, dist_data.shape[1])), columns=dist_data.columns)])
            
            for sat_name, sat_func in self.saturation_models.items():
                result_data = pd.DataFrame({"X": dist_data["X"]})
                for col in dist_data.columns[1:]:
                    mult = dist_data[col].values
                    for i in range(208):
                        mult = np.vstack([mult, self.lagpad(mult, i)])
                    transformed_data = np.nansum(mult * self.data.iloc[:, 1:].values, axis=1)
                    result_data[f"{col}_{sat_name}"] = sat_func(transformed_data)
                
                result_data.iloc[:, 1:] = result_data.iloc[:, 1:].apply(self.ztran, axis=0)
                result_data["Week"] = self.dates
                self.data = result_data
                correlations = self.compute_correlations()
                correlations["Vehicle"] = f"{dist_name}_{sat_name}"
                final_results.append(correlations)
        
        self.results = pd.concat(final_results, ignore_index=True)
        return self

