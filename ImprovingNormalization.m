% Define parameters
N = 1000; % Number of samples
scale = 1; % Scale parameter (fixed)
shape_values = [1, 2, 3]; % Different shape parameters


% Plot combined histograms for each shape parameter
% Define parameters
N = 1000; % Number of samples
scale = 1; % Scale parameter (fixed)
shape = 3; % Shape parameter (fixed)
% Generate normally distributed data for comparison
data = wblrnd(shape, scale, N, 1);
% Initialize vector to store skewness values

% Initialize vector to store skewness values
num_iterations = 8; % Number of iterations
skewness_values = zeros(num_iterations,1);
skewness_values_GP = zeros(num_iterations,5);
Ni = zeros(num_iterations,1);

% Define the range for random indices
lower_bound = 1;
upper_bound = 1000;


% Loop to progressively select smaller samples and calculate skewness
for i = 1:num_iterations
    % Split sample dimension in half
    if i == 1
        random_index = randi([lower_bound, upper_bound], 1, N);
    else
        random_index = randi([lower_bound, upper_bound], 1, round(N / 2));
    end
    
    sample = data(random_index);
    N = length(sample);
    Ni(i) = N;
    
    % Generate normally distributed data for comparison using KCDF_Estim
    [y_kcdf, y_normdata, bgk_band_y, y_Supplement] = KCDF_Estim(sample, 'Tria');
    skewness_values(i) = skewness(y_normdata);
    
    % Generate normally distributed data for comparison using KCDF_EstimGP
    kernel_functions = {'SquaredExponential', 'Matern32', 'Matern52', 'RationalQuadratic', 'Exponential'};
    for j = 1:5
    kernel_function = kernel_functions{j};
    [y_kcdf_GP, y_normdata_GP, bgk_band_y_GP, y_Supplement_GP] = KCDF_EstimGP(sample, 'Tria', kernel_function);
    skewness_values_GP(i,j) = skewness(y_normdata_GP);
    disp(j)
    end
    disp(i)
end

% Define colors
colors = { 'red', 'green', 'blue', 'cyan',
'magenta', 'yellow', 'black'}% Plot the skewness values for both methods
figure;

% Plot each column of skewness_values_GP against N
hold on;
for i = 1:size(skewness_values_GP, 2)
    plot(Ni, skewness_values_GP(:, i),colors{i},'LineWidth', 2);
end
plot(Ni, skewness_values, 'black', 'LineWidth', 2);
xlabel('Sample Size');
ylabel('Skewness');
title('Evolution of Skewness with Decreasing Sample Size');
legend('SquaredExponential', 'Matern32', 'Matern52', 'RationalQuadratic', 'Exponential' ,'Linear Interpolation');
% Plot the skewness values for both methods


figure;
% Plot each column of skewness_values_GP against N
hold on;
for i = 1:size(skewness_values_GP, 2)
    plot(Ni, skewness_values_GP(:, i),'LineWidth', 2, Color=colors{i});
end
plot(Ni, skewness_values, 'b-', 'LineWidth', 2);

xlabel('Sample Size');
ylabel('Skewness');
title('Evolution of Skewness with Decreasing Sample Size');
legend('Linear Interpolation', 'SquaredExponential', 'Matern32', 'Matern52', 'RationalQuadratic', 'Exponential' );

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%second experiment (real data)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
data=target;
N=1000;
% Loop to progressively select smaller samples and calculate skewness
for i = 1:num_iterations
    % Split sample dimension in half
    if i == 1
        random_index = 1:1:N;
    else
        random_index = 1:1:round(N/2);
    end
    
    sample = data(random_index);
    N = length(sample);
    Ni(i) = N;
    
    % Generate normally distributed data for comparison using KCDF_Estim
    [y_kcdf, y_normdata, bgk_band_y, y_Supplement] = KCDF_Estim(sample, 'Tria');
    skewness_values(i) = skewness(y_normdata);
    
    % Generate normally distributed data for comparison using KCDF_EstimGP
    kernel_functions = {'SquaredExponential', 'Matern32', 'Matern52', 'RationalQuadratic', 'Exponential'};
    for j = 1:5
    kernel_function = kernel_functions{j};
    [y_kcdf_GP, y_normdata_GP, bgk_band_y_GP, y_Supplement_GP] = KCDF_EstimGP(sample, 'Tria', kernel_function);
    skewness_values_GP(i,j) = skewness(y_normdata_GP);
    disp(j)
    end
    disp(i)
end

[skewness_values_GP,skewness_values]
colors = { 'red', 'green', 'blue', 'cyan','magenta', 'yellow', 'black'}%
figure;
% Plot each column of skewness_values_GP against N
hold on;
for i = 1:size(skewness_values_GP, 2)
    plot(Ni(1:5), skewness_values_GP(1:5, i), 'LineWidth', 1, 'Color', colors{i});
end
plot(Ni(1:5), skewness_values(1:5), 'black', 'LineWidth', 2);

xlabel('Sample Size');
ylabel('Skewness');
title('Evolution of Skewness with Decreasing Sample Size');
legend('SquaredExponential', 'Matern32', 'Matern52', 'RationalQuadratic', 'Exponential', "Linear Interpolation");




















