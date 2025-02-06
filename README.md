# Adstock Modeling with Carryover and Saturation Effects

## Overview
This repository provides an implementation of an **Adstock Model** that applies **carryover and saturation effects** to media variables in marketing mix modeling. The goal is to estimate the impact of marketing investments over time by utilizing different statistical distributions and saturation functions.

This model allows:
- **Custom distribution selection** for carryover effects.
- **Multiple saturation models** for diminishing returns.
- **Automated correlation computation** to determine the best transformation for each variable.

## Why Use Different Distributions for Carryover and Saturation?
The model applies different probability distributions to model **carryover effects** (how media impact persists over time) and **saturation effects** (how media impact diminishes at higher spend levels).

### Carryover Effects (Distributions)
Each carryover function represents a different way media investments can decay over time:

1. **Weibull Distribution:**
   - Provides a flexible shape for carryover, accommodating **fast and slow decay**.
   - Suitable when media impact **peaks after some delay** before gradually decaying.
   
2. **Lognormal Distribution:**
   - Skewed-right distribution modeling **longer-term carryover**.
   - Useful for campaigns with **a long tail of influence** (e.g., brand-building efforts).
   
3. **Logistic Distribution:**
   - Models media impact that grows **rapidly at first** before stabilizing.
   - Effective for campaigns with **immediate but temporary impact**.
   
4. **Negative Exponential Distribution:**
   - Assumes **instantaneous peak** followed by rapid decay.
   - Suitable for **short-term, high-frequency media** (e.g., social media ads).
   
5. **Chi-Square Distribution:**
   - Models **steep initial impact** followed by a long tail.
   - Useful for **awareness campaigns** that decline slowly over time.

### Saturation Models
Saturation models define how media effectiveness **reduces** as investments increase. We include four types:

1. **Linear Saturation:**
   - No saturation, impact grows **linearly** with spend.
   - Best for **low spend levels where saturation is minimal**.
   
2. **Diminishing Returns:**
   - Modeled as: `impact = spend / (1 + spend)`
   - Captures **gradual saturation**, meaning **diminishing incremental impact**.
   
3. **Logarithmic Saturation:**
   - Modeled as: `impact = log(1 + spend)`
   - Effective for **steep saturation at high spend levels**.
   
4. **Exponential Saturation:**
   - Modeled as: `impact = 1 - exp(-spend)`
   - Represents **rapid initial impact that quickly saturates**.

## How It Works
1. **Load your media data** including marketing investments and sales/equity.
2. **Define carryover distributions** using parameters for each function.
3. **Specify saturation models** to test different diminishing return effects.
4. **Run adstock transformation** for all carryover and saturation combinations.
5. **Compute correlations** between transformed media variables and sales.
6. **Identify the best performing combination** for each variable.

## Installation
To install required dependencies:

For Python:
```bash
pip install -r requirements.txt
```

For R:
```r
install.packages(readLines("requirements.R"))
```

## Example Usage
```python
from adstock_model import AdstockModel

# Define saturation models
saturation_models = {
    "linear": lambda x: x,
    "diminishing": lambda x: x / (1 + x),
    "log": lambda x: np.log1p(x),
    "exponential": lambda x: 1 - np.exp(-x)
}

# Define carryover distributions
distribution_params = {
    "Weibull": {"param1": [0.1, 5, 0.1], "param2": [1, 70, 1]},
    "Lognormal": {"param1": [1, 20, 1], "param2": [10, 20, 1]},
    "Logistic": {"param1": [1, 16, 1], "param2": [1, 16, 1]},
    "NegExp": {"param1": [0.1, 1, 0.1], "param2": [None, None, None]},
    "ChiSquare": {"param1": [1, 70, 1], "param2": [None, None, None]},
}

# Load data
data = pd.read_csv("media_data.csv")
dates = data["Date"]

# Initialize and run the model
adstock_model = AdstockModel(data, dates, saturation_models, distribution_params)
adstock_model.run_adstock()
adstock_model.results.to_csv("adstock_results.csv", index=False)
```

## Output
The model produces a CSV file containing:
- **Transformed media variables** with carryover and saturation applied.
- **Correlation values** to measure effectiveness.
- **Best-performing media transformations**.

## Contribution
Feel free to contribute by submitting pull requests for:
- Adding new carryover functions.
- Enhancing optimization techniques.
- Improving computational performance.

## License
MIT License.


### License
This project is licensed under MIT License. See LICENSE for details.<img width="1240" alt="Scrresnhot" src="https://github.com/user-attachments/assets/8450edd1-2710-4f44-b930-e203b4ac0280">
