%==========================================================================
%	INTRODUCTION
%==========================================================================
%	The main file of the repository is the script "PlotComparisonMF_WMF.m" 
%       The file contains multiple functions referenced below.
%	The script generates a the comparison between the Multifidelity Gaussian Process
%       and the Warped Multifidelity Gaussian process. 

%==========================================================================
%	DATASETS 
%==========================================================================

%       The repository is also equipped with 2 datasets that can be used 
%       to try the models.      

%       1) complete_dataset.csv : It contains 75 time-series of wind-speed data from 2022 
%          from ARPA Lombardia monitoring stations. The data can be also download using the 
           R package ARPALData [5]
  
        2) d_match.csv: It contains the names of the matched time-series and their correlations 
           coefficients, see the paper on clustering experiment. For simplicity only 
           one match is provided.

%==========================================================================
%	OVERVIEW
%==========================================================================

1. This MATLAB script is designed to perform various tasks related to the interpolation of missing sequence of skewed time series. 
-It reads and processes data from CSV files, performs Gaussian Process Regression, and evaluates two different models 
(MFGP and WMFGP) to impute missing values in the time series data. Additionally, 
it calculates errors and generates plots for visualization.

Step-by-Step Explanation:

2. Setting the Working Directory:
   - The script starts by setting the working directory to a specific path.
3. Loading Data:
   - It reads two CSV files: "complete_dataset.csv" and "d_match.csv" into MATLAB tables. These files contain the time series data.
4. Initializing Variables:
   - It defines an array of missing sequence lengths and selects the fifth element as the missing_length.
5. Iterating over Random Time Series (if desired):
   - A for loop iterates over 75 random time series, calculating a discrepancy value for each pair of time series.
6. Gaussian Process Regression:
   - A section of the script performs Gaussian Process Regression on a selected chunk of the time series data, imputing missing values.
7. Model 2 - MFGP:
   - The script defines and optimizes a Multifidelity Gaussian Process model (MFGP) with specific hyperparameters.
   - It predicts missing values using this model.
8. Model 3 - WMFGP:
   - Another section estimates a Kernel Cumulative Distribution Function (KCDF) for low and high fidelity data.
   - It defines and optimizes a Warped Multifidelity Gaussian Process model for WMFGP, making predictions and mapping them back to the original scale.
9. Calculating Mean Absolute Error (MAE) and Root Mean Square Error (RMSE):
   - The script calculates MAE and RMSE for different prediction models, including Gaussian Process Regression, MFGP, and WMFGP.
10. Plotting:
    - Several plots are generated to visualize the observed data, test data, and predictions from different models. 
     The plots include the original data, interpolated data, and predictions,
      along with error bounds.

%==========================================================================
%	NOTE
%==========================================================================

% The "step" parameter controls the sampling frequency of high-fidelity (HF) data points in multifidelity calculations.
% By default, it is set to 3, which means one HF data point is used every 3 observations.
% Users can adjust this parameter based on their specific needs and dataset characteristics.
% If models fail to converge, you can experiment with different "step" values.
% 
% - A value of 1 (step=1) corresponds to using HF data at every observation point. 
%   This dense sampling may be necessary for some datasets.
% 
% - Values between 1 and 5 typically provide good results without a significant loss in performance.
%   For example, step=3 means one HF data point is used every 3 observations.
% 
% - Keep in mind that the choice of "step" should be data-dependent and case-specific.
%   It may require experimentation to find the optimal value that leads to convergence.
% 
% In our research paper, we consistently used a full sampling design (step=1) for HF data 
% to ensure result consistency. However, for practical applications, you can adjust the "step" 
% parameter to balance computation time and accuracy according to your specific dataset.
% 



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Warping is performed using the non-parametric KCDE method presented in Refs. [1], [2] below. 
% The KCDE method is based on the kernel density estimation with the BGK bandwith.
% The repository has been tested under the following versions: Matlab 2023a
%==========================================================================
%   FUNCTIONS USED: Gen_Lookup, kcde, KCDF_Estim, KernBW, kernel_i, Kernel_invNS,
%                   likelihood, sq_dist, minimize, f_H, boundedline
%--------------------------------------------------------------------------
%   The functions: Gen_Lookup, kcde, KCDF_Estim, KernBW, kernel_i, Kernel_invNS 
%   have been written by Andreas Pavlides.
%--------------------------------------------------------------------------
%   The function kde.m has been written by Z. Botev and downloaded from 
%	Mathworks FIle Exchange. 
%-------------------------------------------------------------------------
%   The function: likelihood, sq_dist, minimize, f_H, boundedline have been
%    written by Maziar Raissi  and are available also at the following
%   gitHub repository  https://github.com/maziarraissi/TutorialGP/tree/master/Multifidelity_GP/ 
%   They are publicly available under the following licence:
%MIT License
Copyright (c) 2018 maziarraissi

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
%==========================================================================
%	REFERENCES 
%==========================================================================
% When using this repository please cite the following: 
%
%   1. Agou VD, Pavlides A, Hristopulos DT. Spatial Modeling of 
%   Precipitation Based on Data-Driven Warping of Gaussian Processes. 
%   Entropy. 2022; 24(3):321. https://doi.org/10.3390/e24030321 
%
%   2. Andrew Pavlides, Vasiliki D. Agou, Dionissios T. Hristopulos,
%   Non-parametric kernel-based estimation and simulation of precipitation
%   amount, Journal of Hydrology, Volume 612, Part A, 2022, 127988,
%   https://doi.org/10.1016/j.jhydrol.2022.127988
%
%   3. D. T. Hristopulos, Random Fields for Spatial Data Modeling: 
%   A Primer for Scientists and Engineers. Springer Netherlands, Dordrecht, 2020,
%	ISBN  978-94-024-1918-4,
%	doi: 10.1007/978-94-024-1918-414,
%	URL: https://doi.org/10.1007/978-94-024-1918-4_14.
% 
%   4. For the KDE method with BGK bandwidth cite the paper: 
%      Botev, Z.I., Grotowski, J.F., Kroese, D., 2010. Kernel density estimation via diffusion.
%      Ann. Stat. 38 (5), 2916–2957.
%
%   5. Maranzano, P., Algieri, A., 2023. ARPALData: 
%      an R package for retrieving and analyzing air quality and 
%      weather data from ARPA Lombardia (Italy)

%   6. Colombo, P.,              



