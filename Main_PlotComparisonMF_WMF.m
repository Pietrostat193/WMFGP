% Set working directory
%cd('C:\Users\2692812C\OneDrive - University of Glasgow\Desktop\projects\UsePaoloData\2022\')

% Load data
data = readtable("complete_dataset.csv", 'PreserveVariableNames', true);
data = data(:, 2:end);
match = readtable("d_match.csv");

% Define missing sequence lengths
missing_lengths = [24, 48, 72, 96, 192];
missing_length = missing_lengths(5);





% Initialize iteration and discrepancy array
iter = 1;
dis = zeros(75, 1);

% Iterate over 75 random time series
for i = 1:75
    random_ts_idx = i;
    selected_ts = match(random_ts_idx,:);
    target = table2array(data(:, selected_ts.V2{1}));
    surrogate = table2array(data(:, selected_ts.V3{1}));
    dis(i) = mean(abs(target - surrogate));
end
%This object can be helpfull to understand 
% what type time-series we have and their discrepancy
comparison=[dis,match.V1,(1:75)'];
% Example: High discrepancy between LF and HF and High Correlation
% Random_index can be manually selected or randomized.
% The below index is selected.
random_ts_idx = 8;
selected_ts = match(random_ts_idx,:);
target = table2array(data(:, selected_ts.V2{1}));
surrogate = table2array(data(:, selected_ts.V3{1}));

% Random position of the Gap
seed = 12345 + (iter * 6);
rng(seed);
start_idx = randi([500, length(target) - 400]);
end_idx = start_idx + missing_length - 1;

% Select a chunk of Target and Surrogate
% 200 is the number of data after and before 
% the missing sequence used. It can be change 
% if desired.
a = target((start_idx - 200):(end_idx + 200));
b = surrogate((start_idx - 200):(end_idx + 200));

% Calculate correlation in the chunk of the surrogate
cor = corr(a, b);

% Create a missing sequence
missing_data = target;
test = target(start_idx:end_idx);
missing_data(start_idx:end_idx) = NaN;

% Perform Gaussian Process Regression
Y = missing_data((start_idx - 200):(end_idx + 200));
observed_data = ~isnan(missing_data);
X = (start_idx - 200:end_idx + 200)';
ob = ~isnan(Y);

% Model 1: Train the Gaussian Process model
gprMdl = fitrgp(X(ob), Y(ob));

% Predict the missing values
predicted_data = predict(gprMdl, X(~ob));
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Model 2: MFGP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

hyp = [log([1 1 1 1]) 1 -4 -4];
global ModelInfo
ModelInfo.hyp = hyp;
ModelInfo.jitter = eps;

% Define input data for Model 2
X_H = (X(ob) - (start_idx - 201));
ends = length(X_H);
step = 3;

ModelInfo.X_H = X_H(1:step:ends);
ModelInfo.X_L = X - (start_idx - 201);

% Output data for Model 2
ModelInfo.y_L = surrogate(X);
y_H = Y(ob);
ModelInfo.y_H = y_H(1:step:ends);

% Optimization for Model 2
cd 'C:\Users\2692812C\OneDrive - University of Glasgow\Desktop\PietroPhD\5.WMF\NewMethod\Test1D_GPR_public\FolderToShare'
options = optimoptions('fminunc', 'GradObj', 'on', 'Display', 'iter', ...
    'Algorithm', 'trust-region', 'Diagnostics', 'on', 'DerivativeCheck', 'on', ...
    'FinDiffType', 'central');

[hyp, fval, ~, ~, ~, ~] = fminunc(@likelihood, hyp, options);

ModelInfo.hyp = hyp;
x_star = ModelInfo.X_L;
% Prediction
[MF_predicted, var_MF] = predictor_f_H(x_star);

      
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %Model 3: WMFGP
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        kernel = 'Tria';
        % KCDF Estimation for Low Fidelity Data
        [y_L_kcdf, y_L_normdata, bgk_band_y_L, y_L_Supplement] = KCDF_Estim(surrogate, 'Tria');
        [lookup_y_L] = Gen_Lookup(surrogate, bgk_band_y_L, kernel);
        ny_L = y_L_normdata((start_idx-200):(end_idx+200));
        
        % High Fidelity Normalization
        [y_H_kcdf, y_H_normdata, bgk_band_y_H, y_H_Supplement] = KCDF_Estim(target, 'Tria');
        [lookup_y_H] = Gen_Lookup(target, bgk_band_y_H, kernel);
        %Y=missing_data((start_idx-200):(end_idx+200));
        ny_H=y_H_normdata((start_idx-200):(end_idx+200));
        ny_H=ny_H(ob);
        ny_H=ny_H(1:step:ends);
        
       
            hyp = [log([1 1 1 1]) 1 -4 -4];
            ModelInfo.y_L=ny_L;
            %X_H=(X(ob)-(start_idx-201));
            %ModelInfo.X_H=X_H;
            ModelInfo.y_H=ny_H;
            [hyp,fval, ~, ~, ~, ~] = fminunc(@likelihood, hyp, options);
            ModelInfo.hyp=hyp;
            x_star=ModelInfo.X_L;
            %% Predicition
            [WMF_predicted, WMF_var]=predictor_f_H(x_star);
            %   Calculate the quantiles of the normal scores
            znorm_low = WMF_predicted -2 * sqrt(WMF_var);
            znorm_high =WMF_predicted + 2 * sqrt(WMF_var);
            %Backmapping
            [Zhat_warp_HF] = Kernel_invNS(WMF_predicted, lookup_y_H);
            [lower_warps] = Kernel_invNS(znorm_low, lookup_y_H);
            [upper_warps] = Kernel_invNS(znorm_high, lookup_y_H);
            
            

 % Calculate MAE and RMSE
actual_missing_values = target(start_idx:end_idx);
 
mae = mean(abs(predicted_data - actual_missing_values));
rmse = sqrt(mean((predicted_data - actual_missing_values).^2));
        
mae_MF = mean(abs(MF_predicted(~ob) - actual_missing_values));
rmse_MF = sqrt(mean((MF_predicted(~ob) - actual_missing_values).^2));
        
mae_WMF = mean(abs(Zhat_warp_HF(~ob) - actual_missing_values));
rmse_WMF = sqrt(mean((Zhat_warp_HF(~ob) - actual_missing_values).^2));

[mae, mae_MF, mae_WMF]


%% Plot
X_L=min(X(ob)):max(X(ob));
y_L=surrogate(X);
X_H=X_H+(start_idx - 201);
% Plot original data and predictions for Model 1

figure;
subplot(1,2,1)
hold on;
% Plot observed LF data with blue circles
plot(X_L, y_L, 'bo', 'DisplayName', 'Observed LF Data', 'MarkerSize', 3, 'LineWidth', 1.5);
% Plot observed HF data with red squares
plot(X_H(1:step:ends), y_H(1:step:ends), 'rs', 'DisplayName', 'Observed HF Data', 'MarkerSize', 3, 'LineWidth', 1.5);
% Shaded interpolation AREA
x_shaded = [start_idx, end_idx, end_idx, start_idx];
y_shaded = [min(y_L - 10), min(y_L - 10), max(y_L + 10), max(y_L + 10)];
fill(x_shaded, y_shaded, 'y', 'FaceAlpha', 0.2, 'DisplayName', 'Interpolation Area');
legend('Observed LF Data','Observed HF Data', 'Interpolation Area')
ylim([min(y_H - 2), max(y_H + 2)]); % Set
title('Training Data and Interpolation Area');
xlabel('Time-Index');
ylabel('m/s');
grid on; % Show grid lines


subplot(1,2,2)
hold on;
% Plot observed LF data with blue circles
plot(X_L, y_L, 'bo', 'DisplayName', 'Observed LF Data', 'MarkerSize', 3, 'LineWidth', 1.5);

% Plot observed HF data with red squares
plot(X_H(1:step:ends), y_H(1:step:ends), 'rs', 'DisplayName', 'Observed HF Data', 'MarkerSize', 3, 'LineWidth', 1.5);

% Plot missing data as red dots
plot(start_idx:end_idx, actual_missing_values, 'r.', 'DisplayName', 'Test Data', 'MarkerSize', 10);

% Plot MF predictions with cyan dotted line
plot(X_L, MF_predicted, 'c:', 'DisplayName', 'Predicted Data MF', 'LineWidth', 1.5);

% Plot WMF predictions with magenta dash-dot line
%plot(X_L, Zhat_warp_HF, 'm.-', 'DisplayName', 'Predicted Data WMF', 'LineWidth', 1.2);
[l,h(4)] = boundedline(X_L, Zhat_warp_HF, (upper_warps-Zhat_warp_HF)' ,':', 'alpha','cmap','magenta');
% Shaded interpolation AREA
x_shaded = [start_idx, end_idx, end_idx, start_idx];
y_shaded = [min(y_L - 10), min(y_L - 10), max(y_L + 10), max(y_L + 10)];
fill(x_shaded, y_shaded, 'y', 'FaceAlpha', 0.2, 'DisplayName', 'Interpolation Area');

% Adjust axis limits and other properties as needed
xlim([start_idx - 30, end_idx + 40]); % Set the x-axis limits
ylim([min(y_H - 3), max(y_H + 3)]); % Set the y-axis limits
grid on; % Show grid lines

title('Model Predictions');
xlabel('Time-Index');
ylabel('m/s');
legend('Observed LF Data','Observed HF Data', 'Test Data', 'Predicted Data MF', 'Predicted Data WMF', "CI Predicted Data MF",'Interpolation Area');

%%%%%%%%%%%%%%%%%%
% Third Plot for zoom
%%%%%%%%%%%%%%%%%%%%%%
figure;
hold on;
% Plot observed LF data with blue circles
plot(X_L, y_L, 'bo', 'DisplayName', 'Observed LF Data', 'MarkerSize', 3, 'LineWidth', 1.5);

% Plot observed HF data with red squares
plot(X_H(1:step:ends), y_H(1:step:ends), 'rs', 'DisplayName', 'Observed HF Data', 'MarkerSize', 3, 'LineWidth', 1.5);

% Plot missing data as red dots
plot(start_idx:end_idx, actual_missing_values, 'r.', 'DisplayName', 'Test Data', 'MarkerSize', 10);

% Plot MF predictions with cyan dotted line
plot(X_L, MF_predicted, 'c:', 'DisplayName', 'Predicted Data MF', 'LineWidth', 1.5);

% Plot WMF predictions with magenta dash-dot line
%plot(X_L, Zhat_warp_HF, 'm.-', 'DisplayName', 'Predicted Data WMF', 'LineWidth', 1.2);
[l,h(4)] = boundedline(X_L, Zhat_warp_HF, (upper_warps-Zhat_warp_HF)' ,':', 'alpha','cmap','magenta');
% Shaded interpolation AREA
x_shaded = [start_idx, end_idx, end_idx, start_idx];
y_shaded = [min(y_L - 10), min(y_L - 10), max(y_L + 10), max(y_L + 10)];
fill(x_shaded, y_shaded, 'y', 'FaceAlpha', 0.2, 'DisplayName', 'Interpolation Area');

% Adjust axis limits and other properties as needed
xlim([start_idx+50, end_idx-50]); % Set the x-axis limits
ylim([min(y_H - 2), max(y_H + 2)]); % Set the y-axis limits
grid on; % Show grid lines

title('Model Predictions');
xlabel('Time-Index');
ylabel('m/s');
legend('Observed LF Data','Observed HF Data', 'Test Data', 'Predicted Data MF', 'Predicted Data WMF', "CI Predicted Data MF",'Interpolation Area');
