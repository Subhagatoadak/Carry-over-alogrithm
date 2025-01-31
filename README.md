# Adstock Carryover Algorithm for Media Impact Analysis
### Overview
This tool implements an Adstock Carryover Algorithm to measure the short-term and long-term impact of media vehicles on key business metrics (e.g., sales). It enables users to select the best-correlated media vehicle based on probability distributions that capture media retention effects over time.

### Features
##### Adstock Function Selection: Users can choose different Adstock functions (e.g., logistic, exponential) to model media decay.
##### Probability Distribution Modeling: Implements statistical distributions (e.g., logistic distribution) to estimate retention effects.
##### Correlation Evaluation: Automatically evaluates correlations between media impressions, Adstock transformation, and business KPIs.
##### Visualization Dashboard: Provides interactive plots to compare raw media impressions, Adstock-transformed values, and business outcomes.
##### Customizable Parameters: Users can adjust location and scale parameters to optimize media decay behavior.

### How It Works
Select a Media Vehicle – The tool allows users to analyze different media vehicles (e.g., TV, display ads, search ads).
Apply an Adstock Function – Choose an appropriate decay function to transform media impressions into long-term impact values.
Configure Distribution Parameters – Adjust retention probabilities using location and scale parameters.
Analyze Correlation Metrics – Evaluate correlation strength between transformed media values and business KPIs.
Interpret Visualizations – Leverage the dashboard to understand short- and long-term media carryover effects.

Modify Adstock parameters in the UI to analyze different retention behaviors.
Compare correlation results to identify the most impactful media vehicle.
Export results for further modeling in Marketing Mix Models (MMM) or other analytics workflows.
Example Output
The dashboard provides insights via:

Correlation metrics (e.g., -14% with variable, -16% with Adstock).
Time series plots showing media impressions vs. transformed Adstock values.
Distribution plots to visualize media decay effects.
Contributing
Contributions are welcome! Please submit a pull request or open an issue to suggest improvements.

### License
This project is licensed under MIT License. See LICENSE for details.<img width="1240" alt="Scrresnhot" src="https://github.com/user-attachments/assets/8450edd1-2710-4f44-b930-e203b4ac0280">
