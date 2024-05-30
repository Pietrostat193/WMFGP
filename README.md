# Warped Multifidelity Gaussian Process Models for Environmental Data Imputation

## Introduction
The main file of this repository is the script `PlotComparisonMF_WMF.m`. This script generates a comparison between the Multifidelity Gaussian Process (MFGP) and the Warped Multifidelity Gaussian Process (WMFGP).

## Datasets
The repository includes two datasets:

1. `complete_dataset.csv`: Contains 75 time-series of wind-speed data from 2022 from ARPA Lombardia monitoring stations. The data can also be downloaded using the R package ARPALData.
2. `d_match.csv`: Contains the names of the matched time-series and their correlation coefficients. Only one match is provided for simplicity.

## Overview

This MATLAB script is designed to perform various tasks related to the interpolation of missing sequences in skewed time series. It reads and processes data from CSV files, performs Gaussian Process Regression, and evaluates two different models (MFGP and WMFGP) to impute missing values in the time series data. Additionally, it calculates errors and generates plots for visualization.

### Step-by-Step Explanation

1. **Setting the Working Directory:**
   - The script starts by setting the working directory to a specific path.
2. **Loading Data:**
   - It reads two CSV files: `complete_dataset.csv` and `d_match.csv` into MATLAB tables. These files contain the time series data.
3. **Initializing Variables:**
   - It defines an array of missing sequence lengths and selects the fifth element as the missing_length.
4. **Iterating over Random Time Series (if desired):**
   - A for loop iterates over 75 random time series, calculating a discrepancy value for each pair of time series.
5. **Gaussian Process Regression:**
   - Performs Gaussian Process Regression on a selected chunk of the time series data, imputing missing values.
6. **Model 2 - MFGP:**
   - Defines and optimizes a Multifidelity Gaussian Process model (MFGP) with specific hyperparameters and predicts missing values using this model.
7. **Model 3 - WMFGP:**
   - Estimates a Kernel Cumulative Distribution Function (KCDF) for low and high fidelity data, defines and optimizes a Warped Multifidelity Gaussian Process model (WMFGP), making predictions and mapping them back to the original scale.
8. **Calculating Mean Absolute Error (MAE) and Root Mean Square Error (RMSE):**
   - Calculates MAE and RMSE for different prediction models, including Gaussian Process Regression, MFGP, and WMFGP.
9. **Plotting:**
    - Generates several plots to visualize the observed data, test data, and predictions from different models. The plots include the original data, interpolated data, and predictions, along with error bounds.

## Note

- The `step` parameter controls the sampling frequency of high-fidelity (HF) data points in multifidelity calculations. By default, it is set to 3, which means one HF data point is used every 3 observations.
- Users can adjust this parameter based on their specific needs and dataset characteristics. Experimentation may be necessary to find the optimal value that leads to convergence.
- In the research paper, a full sampling design (step=1) was used for HF data to ensure result consistency. For practical applications, you can adjust the `step` parameter to balance computation time and accuracy according to your specific dataset.

## Functions Used

- **Custom Functions by Andreas Pavlides:**
  - `Gen_Lookup`, `kcde`, `KCDF_Estim`, `KernBW`, `kernel_i`, `Kernel_invNS`
- **From Mathworks File Exchange:**
  - `kde.m` by Z. Botev
- **From Maziar Raissi's GitHub Repository:**
  - `likelihood`, `sq_dist`, `minimize`, `f_H`, `boundedline` (MIT License)

## References

When using this repository, please cite the following references:

1. Agou VD, Pavlides A, Hristopulos DT. "Spatial Modeling of Precipitation Based on Data-Driven Warping of Gaussian Processes." Entropy. 2022; 24(3):321. [https://doi.org/10.3390/e24030321](https://doi.org/10.3390/e24030321)
2. Pavlides A, Agou VD, Hristopulos DT. "Non-parametric kernel-based estimation and simulation of precipitation amount." Journal of Hydrology, 2022, [https://doi.org/10.1016/j.jhydrol.2022.127988](https://doi.org/10.1016/j.jhydrol.2022.127988)
3. Hristopulos DT. "Random Fields for Spatial Data Modeling: A Primer for Scientists and Engineers." Springer, 2020, [https://doi.org/10.1007/978-94-024-1918-4_14](https://doi.org/10.1007/978-94-024-1918-4_14)
4. Botev, Z.I., Grotowski, J.F., Kroese, D. "Kernel density estimation via diffusion." Ann. Stat. 38 (5), 2010, 2916–2957.
5. Maranzano, P., Algieri, A. "ARPALData: an R package for retrieving and analyzing air quality and weather data from ARPA Lombardia (Italy)", 2023.
