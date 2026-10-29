function economic_analysis_fourgp_gp_glbl_CA_second()

% Load data
data = readtable('Economic_data.xlsx');

% Display column headers to check valid variable names
disp(data.Properties.VariableNames);

% Assuming we have identified the correct variable names from the display
% Replace these with actual column names as displayed by the above disp function
gdp_percapita_var = 'ln_gdp_percap';  
labour_force_var = 'ln_labour_force';          
capital_var = 'ln_physical_capital';               
human_capital_var = 'ln_human_capital';   
RD_stock_var = 'ln_rnd_domestic';            
RD_spillovers_var = 'ln_rnd_imp_spill';  
infrastructure_var = 'ln_infrastructure';
urbanization_var = 'ln_urbanisation'; 
region_var = 'region'; 
year_var = 'year';                    
country_var = 'country';

% Combine "Europe" and "North America and Oceania" into one region
data.(region_var) = regexprep(data.(region_var), 'North America and Oceania|Europe', 'Europe and North America and Oceania');

% Define the four regions of interest
regions_of_interest = {'Africa', 'Asia', 'Europe and North America and Oceania', 'Latin America and the Caribbean'};

% Filter data to include only the specified regions
data = data(ismember(data.(region_var), regions_of_interest), :);

% Compute global cross-sectional average for GDP per capita
global_avg_gdp_percapita = nanmean(data.(gdp_percapita_var));

% Group data by region
regions = unique(data.(region_var));
numRegions = length(regions);

% Initialize result storage
results = struct();
cce_estimates = cell(10, numRegions); % To store CCE estimates for all regions
ccemg_estimates = nan(10, numRegions); % To store CCEMG estimates for each region
gdp_percapita_data = struct(); % To store GDP per capita data for each region

for r = 1:numRegions
    regionName = regions{r};
    sanitizedRegionName = matlab.lang.makeValidName(regionName);
    
    regionData = data(strcmp(data.(region_var), regionName), :);
    
    % Get unique years and countries for reshaping
    years = unique(regionData.(year_var));
    countries = unique(regionData.(country_var));
    T = length(years);
    N = length(countries);
    
    % Initialize matrices
    Y = nan(T, N);
    X = nan(T, N, 11); % Increased number of covariates to 11
    
    for i = 1:N
        countryData = regionData(strcmp(regionData.(country_var), countries{i}), :);
        
        for t = 1:T
            yearData = countryData(countryData.(year_var) == years(t), :);
            if ~isempty(yearData)
                Y(t, i) = yearData.(gdp_percapita_var);
                X(t, i, 1) = yearData.(labour_force_var);
                X(t, i, 2) = yearData.(capital_var);
                X(t, i, 3) = yearData.(human_capital_var);
                X(t, i, 4) = yearData.(RD_stock_var);
                X(t, i, 5) = yearData.(RD_spillovers_var);
                X(t, i, 6) = yearData.(infrastructure_var);
                X(t, i, 7) = yearData.(urbanization_var);
                X(t, i, 8) = global_avg_gdp_percapita; % Add global average GDP per capita as a covariate
            end
        end
    end
    
    % Compute cross-sectional averages for additional covariates
    cross_sectional_averages = nanmean(Y, 2);
    cross_sectional_averages_lf = nanmean(squeeze(X(:, :, 1)), 2);
    cross_sectional_averages_pc = nanmean(squeeze(X(:, :, 2)), 2);
    
    % Add cross-sectional averages to X matrix as new covariates
    for t = 1:T
        X(t, :, 9) = cross_sectional_averages(t);
        X(t, :, 10) = cross_sectional_averages_lf(t);
        X(t, :, 11) = cross_sectional_averages_pc(t);
    end
    
    % Perform CCE and CCEMG estimation
    [cce_params, ccemg_params, ccemg_sd, ccemg_sdpes, factors, var_ccemg, varpes_ccemg] = CCE_estimator(Y, X, T, N);
    
    % Store results
    results.(sanitizedRegionName).cce_params = cce_params;
    results.(sanitizedRegionName).ccemg_params = ccemg_params;
    results.(sanitizedRegionName).ccemg_sd = ccemg_sd/(sqrt(N*T));
    results.(sanitizedRegionName).ccemg_sdpes = ccemg_sdpes/(sqrt(N));
    results.(sanitizedRegionName).factors = factors;
    results.(sanitizedRegionName).var_ccemg = var_ccemg/(N*T);
    results.(sanitizedRegionName).varpes_ccemg = varpes_ccemg/N;
    
    % Calculate 95% confidence intervals for CCEMG parameters using both methods
    ci_lower = ccemg_params - 1.96 * ccemg_sd/(sqrt(N*T));
    ci_upper = ccemg_params + 1.96 * ccemg_sd/(sqrt(N*T));
    results.(sanitizedRegionName).ccemg_ci = [ci_lower, ci_upper];
    
    ci_lower_pes = ccemg_params - 1.96 * ccemg_sdpes/(sqrt(N));
    ci_upper_pes = ccemg_params + 1.96 * ccemg_sdpes/(sqrt(N));
    results.(sanitizedRegionName).ccemg_ci_pes = [ci_lower_pes, ci_upper_pes];
    
    % Conduct hypothesis tests for CCEMG parameters
    t_stats = ccemg_params ./ (ccemg_sd/(sqrt(N*T)));
    t_stats_pes = ccemg_params ./ (ccemg_sdpes/(sqrt(N)));
    p_values = 2 * (1 - tcdf(abs(t_stats), N - 1)); % Two-tailed test
    p_values_pes = 2 * (1 - tcdf(abs(t_stats_pes), N - 1)); % Two-tailed test
    critical_value = tinv(0.975, N - 1); % 95% confidence level

    results.(sanitizedRegionName).t_stats = t_stats;
    results.(sanitizedRegionName).t_stats_pes = t_stats_pes;
    results.(sanitizedRegionName).p_values = p_values;
    results.(sanitizedRegionName).p_values_pes = p_values_pes;
    results.(sanitizedRegionName).critical_value = critical_value;
    
    % Analyze growth patterns
    gdp_growth_rate = diff(log(Y)); % Log difference to compute growth rate
    mean_growth_rate = mean(gdp_growth_rate, 'omitnan');
    std_growth_rate = std(gdp_growth_rate, 'omitnan');
    
    results.(sanitizedRegionName).gdp_growth_rate = gdp_growth_rate;
    results.(sanitizedRegionName).mean_growth_rate = mean_growth_rate;
    results.(sanitizedRegionName).std_growth_rate = std_growth_rate;
    
    % Compute summary statistics for CCE and CCEMG
    results.(sanitizedRegionName).cce_mean = mean(cce_params, 2, 'omitnan');
    results.(sanitizedRegionName).cce_std = std(cce_params, 0, 2, 'omitnan');
    
    % Store CCE and CCEMG estimates for all regions
    for v = 1:11
        cce_estimates{v, r} = cce_params(v, :);
        ccemg_estimates(v, r) = ccemg_params(v);
    end
    
    % Store GDP per capita data for this region
    gdp_percapita_data.(sanitizedRegionName) = nanmean(Y, 2);

    % Variable names for plots
    variable_names = {labour_force_var, capital_var, human_capital_var, ...
                      RD_stock_var, RD_spillovers_var, infrastructure_var, ...
                      urbanization_var, 'Global Average GDP', ...
                      'Cross-sectional GDP', 'Cross-sectional Labour Force', ...
                      'Cross-sectional Human Capital'};
end

% Determine y-axis limits for each variable across all regions
y_limits = zeros(11, 2);
for v = 1:11
    all_estimates = cell2mat(cce_estimates(v, :));
    y_limits(v, :) = [min(all_estimates, [], 'all'), max(all_estimates, [], 'all')];
end

% Plot scatter plots of CCE estimates for each region in one figure for each variable
for v = 1:11
    figure;
    for r = 1:numRegions
        subplot(ceil(numRegions/2), 2, r);
        regionName = regions{r};
        scatter(1:length(cce_estimates{v, r}), cce_estimates{v, r}, 'filled');
        hold on;
        % Add dashed line for CCEMG estimate
        yline(ccemg_estimates(v, r), '--r', 'LineWidth', 2, 'DisplayName', 'CCEMG Estimate');
        title(['Region: ', regionName]);
        xlabel('Country');
        ylabel(['CCE Estimate of ', variable_names{v}]);
        ylim(y_limits(v, :)); % Set uniform y-axis limits
        legend('show');
        grid on;
        hold off;
    end
    sgtitle(['CCE Estimates for ', variable_names{v}, ' in All Regions']);
end

% Plot GDP per capita over time for all regions in the same plot
figure;
hold on;
for r = 1:numRegions
    regionName = regions{r};
    sanitizedRegionName = matlab.lang.makeValidName(regionName);
    plot(years, gdp_percapita_data.(sanitizedRegionName), 'LineWidth', 2, 'DisplayName', regionName);
end
title('GDP per Capita over Time for All Regions');
xlabel('Year');
ylabel('GDP per Capita');
legend('show');
grid on;
hold off;

% Display results summary
for r = 1:numRegions
    regionName = regions{r};
    sanitizedRegionName = matlab.lang.makeValidName(regionName);
    
    fprintf('Results for %s:\n', regionName);
    fprintf('CCE Mean Estimates:\n');
    disp(results.(sanitizedRegionName).cce_mean);
    fprintf('CCE Standard Deviations:\n');
    disp(results.(sanitizedRegionName).cce_std);
    fprintf('CCEMG Mean Estimates:\n');
    disp(results.(sanitizedRegionName).ccemg_params);
    fprintf('CCEMG Standard Deviations (using Pesaran method):\n');
    disp(results.(sanitizedRegionName).ccemg_sd);
    fprintf('CCEMG Standard Deviations (using alternative method):\n');
    disp(results.(sanitizedRegionName).ccemg_sdpes);
    fprintf('Variance of CCEMG Estimates (using Pesaran method):\n');
    disp(results.(sanitizedRegionName).var_ccemg);
    fprintf('Variance of CCEMG Estimates (using alternative method):\n');
    disp(results.(sanitizedRegionName).varpes_ccemg);
    fprintf('95%% Confidence Intervals for CCEMG Estimates:\n');
    disp(results.(sanitizedRegionName).ccemg_ci);
    fprintf('95%% Confidence Intervals for CCEMG Estimates (Pesaran method):\n');
    disp(results.(sanitizedRegionName).ccemg_ci_pes);
    fprintf('T-Statistics for CCEMG Estimates (Proposed method):\n');
    disp(results.(sanitizedRegionName).t_stats);
    fprintf('T-Statistics for CCEMG Estimates (Pesaran method):\n');
    disp(results.(sanitizedRegionName).t_stats_pes);
    fprintf('P-Values for CCEMG Estimates (Proposed method):\n');
    disp(results.(sanitizedRegionName).p_values);
    fprintf('P-Values for CCEMG Estimates (Pesaran method):\n');
    disp(results.(sanitizedRegionName).p_values_pes);
    fprintf('Critical Value for t-test: %f\n', results.(sanitizedRegionName).critical_value);
    fprintf('Mean GDP Growth Rate: %f\n', results.(sanitizedRegionName).mean_growth_rate);
    fprintf('Standard Deviation of GDP Growth Rate: %f\n', results.(sanitizedRegionName).std_growth_rate);
    fprintf('\n');
end

% Analyze and compare the results obtained from Pesaran's method and the proposed method
compare_results(results, regions, variable_names);

% Plot scatter plot of CCE estimates for "Europe and North America and Oceania" region and "ln_rnd_domestic" variable
plot_scatter_europe_ln_rnd_domestic(results, cce_estimates, ccemg_estimates, regions, 'Europe and North America and Oceania', RD_stock_var);

end

% Function for CCE estimation
function [cce_params, ccemg_params, ccemg_sd, ccemg_sdpes, factors, var_ccemg, varpes_ccemg] = CCE_estimator(Y, X, T, N)
    % Assuming we use cross-sectional means as factors
    factors = mean(X, 2);
    factors = reshape(factors, T, []);
    
    % Initialize beta_i to store coefficients per cross-section
    num_regressors = size(X, 3);
    beta_i = nan(num_regressors, N);  % Coefficients for each regressor
    
    for i = 1:N
        X_i = squeeze(X(:, i, :));
        Y_i = Y(:, i);
        
        % Check if X_i is rank deficient and handle it
        if rank(X_i) < num_regressors
            warning('X_i is rank deficient for country %d.', i);
            continue;
        end
        
        % Fit model without intercept if model simplicity is desired
        model_i = fitlm(X_i, Y_i, 'Intercept', false);
        if size(model_i.Coefficients.Estimate, 1) == num_regressors
            beta_i(:, i) = model_i.Coefficients.Estimate;  % Capture only beta estimates
        else
            warning('Unexpected coefficient size for country %d.', i);
        end
    end
    
    cce_params = beta_i;
    ccemg_params = nanmean(beta_i, 2);  % Average across all cross-sections for CCEMG estimates
    
    % Variance estimation for CCEMG
    covMG1 = zeros(num_regressors, num_regressors); % Covariance matrix (Pesaran)
    covMG4 = zeros(num_regressors, num_regressors); % Covariance for the fourth term
    
    for i = 1:N
        if isnan(beta_i(1, i))
            continue;
        end
        X_i = squeeze(X(:, i, :));
        Y_i = Y(:, i);
        mat1i = beta_i(:,i) - ccemg_params;
        covMG1 = covMG1 + (mat1i)*(mat1i)'; % Estimating covariance 1
        
        e_i = Y_i - X_i * beta_i(:, i);
        x_i = X_i;
        
        m = floor(T^(1/4));
        sumi = 0;
        for t = 1:m
            h_it = e_i(t)' * x_i(t, :)';
            sumit = sumi + h_it * h_it';
            sumi = sumit;
        end
        
        sigmai0 = (1/T)*sumi;
        suminew = 0;
        for j = 1:m
            weightedi = 0;
            for t = j+1:m
                h_it = e_i(t)' * x_i(t, :)';
                h_itj = e_i(t-j)' * x_i(t-j, :)';
                sumitj = suminew + (h_it * h_itj' + h_itj * h_it');
            end
            weighteditj = (1 - (j / (m + 1))) * (1 / T) * sumitj;
            weightedi = weightedi + weighteditj;
        end
    
        Siepsilon = sigmai0 + weightedi; % Estimating covariance 3
        covMG4 = covMG4 + inv((X_i' * X_i) / T) * Siepsilon * inv((X_i' * X_i) / T);
    end

    estcovMG = (T / N) * covMG1 + (1 / N) * covMG4; % Estimated covariance matrix
    estcovMGpes = (1 / N) * covMG1;
    var_ccemg = diag(estcovMG); % Variance of the CCEMG estimator
    varpes_ccemg = diag(estcovMGpes);
    ccemg_sd = sqrt(var_ccemg);    % Standard deviation across all cross-sections for CCEMG estimates
    ccemg_sdpes = sqrt(varpes_ccemg);
end

function plot_scatter_europe_ln_rnd_domestic(results, cce_estimates, ccemg_estimates, regions, target_region, variable_name)
    % Find the index for the target region
    region_idx = find(strcmp(regions, target_region));
    if isempty(region_idx)
        error('Region "%s" not found.', target_region);
    end
    
    % Extract the necessary data
    sanitizedRegionName = matlab.lang.makeValidName(target_region);
    cce_data = cce_estimates{4, region_idx}; % 4 corresponds to 'ln_rnd_domestic'
    ccemg_estimate = ccemg_estimates(4, region_idx);
    ci_proposed = results.(sanitizedRegionName).ccemg_ci(4, :);
    ci_pesaran = results.(sanitizedRegionName).ccemg_ci_pes(4, :);
    
    % Plot the scatter plot
    figure;
    scatter(1:length(cce_data), cce_data, 'filled');
    hold on;
    yline(ccemg_estimate, '--r', 'LineWidth', 2, 'DisplayName', 'CCEMG Estimate');
    yline(ci_proposed(1), ':b', 'LineWidth', 2, 'DisplayName', '95% CI Proposed Lower');
    yline(ci_proposed(2), ':b', 'LineWidth', 2, 'DisplayName', '95% CI Proposed Upper');
    yline(ci_pesaran(1), ':g', 'LineWidth', 2, 'DisplayName', '95% CI Pesaran Lower');
    yline(ci_pesaran(2), ':g', 'LineWidth', 2, 'DisplayName', '95% CI Pesaran Upper');
    title(['CCE Estimates for ', variable_name, ' in ', target_region]);
    xlabel('Country');
    ylabel(['CCE Estimate of ', variable_name]);
    legend('show');
    grid on;
    hold off;
end

function compare_results(results, regions, variable_names)
    numRegions = length(regions);

    % Compare standard deviations and variances
    fprintf('Comparison of Standard Deviations and Variances:\n');
    for v = 1:10
        fprintf('Variable: %s\n', variable_names{v});
        for r = 1:numRegions
            regionName = regions{r};
            sanitizedRegionName = matlab.lang.makeValidName(regionName);
            
            ccemg_sd = results.(sanitizedRegionName).ccemg_sd(v);
            ccemg_sdpes = results.(sanitizedRegionName).ccemg_sdpes(v);
            var_ccemg = results.(sanitizedRegionName).var_ccemg(v);
            varpes_ccemg = results.(sanitizedRegionName).varpes_ccemg(v);
            
            fprintf('Region: %s\n', regionName);
            fprintf('Standard Deviation (Pesaran): %f, Standard Deviation (Proposed): %f\n', ccemg_sd, ccemg_sdpes);
            fprintf('Variance (Pesaran): %f, Variance (Proposed): %f\n', var_ccemg, varpes_ccemg);
        end
        fprintf('\n');
    end

    % Compare confidence intervals
    fprintf('Comparison of 95%% Confidence Intervals:\n');
    for v = 1:10
        fprintf('Variable: %s\n', variable_names{v});
        for r = 1:numRegions
            regionName = regions{r};
            sanitizedRegionName = matlab.lang.makeValidName(regionName);
            
            ci = results.(sanitizedRegionName).ccemg_ci(v, :);
            ci_pes = results.(sanitizedRegionName).ccemg_ci_pes(v, :);
            
            fprintf('Region: %s\n', regionName);
            fprintf('Confidence Interval (Pesaran): [%f, %f]\n', ci(1), ci(2));
            fprintf('Confidence Interval (Proposed): [%f, %f]\n', ci_pes(1), ci_pes(2));
        end
        fprintf('\n');
    end

    % Compare t-statistics
    fprintf('Comparison of T-Statistics:\n');
    for v = 1:10
        fprintf('Variable: %s\n', variable_names{v});
        for r = 1:numRegions
            regionName = regions{r};
            sanitizedRegionName = matlab.lang.makeValidName(regionName);
            
            t_stats = results.(sanitizedRegionName).t_stats(v);
            t_stats_pes = results.(sanitizedRegionName).t_stats_pes(v);
            
            fprintf('Region: %s\n', regionName);
            fprintf('T-Statistic (Proposed): %f, T-Statistic (Pesaran): %f\n', t_stats, t_stats_pes);
        end
        fprintf('\n');
    end
end
