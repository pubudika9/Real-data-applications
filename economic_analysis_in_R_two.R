library(readxl)
library(dplyr)
library(ggplot2)
library(plm)
library(car) # For VIF calculation

aus_wheat_yield <- function() {
  # Load data
  data <- read_excel('All_states_regressors_SPV.csv')
  
  # Display column headers to check valid variable names
  print(colnames(data))
  
  # Assuming we have identified the correct variable names from the display
  # gdp_percapita_var <- 'ln_gdp_percap'
  # labour_force_var <- 'ln_labour_force'
  # capital_var <- 'ln_physical_capital'
  # human_capital_var <- 'ln_human_capital'
  # RD_stock_var <- 'ln_rnd_domestic'
  # RD_spillovers_var <- 'ln_rnd_imp_spill'
  # infrastructure_var <- 'ln_infrastructure'
  # urbanization_var <- 'ln_urbanisation'
  # region_var <- 'region'
  # year_var <- 'year'
  # country_var <- 'country'
  # 
  # Define variables
  lnyield <-  'ln_yields' # log(wheat yield in ton per hectare)
  statistical_area <-  'Region' # Cross-sections
  year <-  'Year' # Time points
  state_var <-  'State'
  
  # Annual cumulative temperature and rainfall variables
  fds1 <-  'FD_s1'
  fds2 <-  'FD_s2'
  fds3 <-  'FD_s3'
  ddlows1 <-  'DD1_s1'
  ddlows2 <-  'DD1_s2'
  ddlows3 <-  'DD1_s3'
  ddmods1 <-  'DD2_s1'
  ddmods2 <-  'DD2_s2'
  ddmods3 <-  'DD2_s3'
  ddhighs1 <-  'DD3_s1'
  ddhighs2 <-  'DD3_s2'
  ddhighs3 <- 'DD3_s3'
  rs0 <-  'r_s0'
  rs1 <-  'r_s1'
  rs2 <-  'r_s2' 
  rs3 <- 'r_s3'
  # Combine "Europe" and "North America and Oceania" into one region
 # data[[state_var]] <- gsub('North America and Oceania|Europe', 'Europe and North America and Oceania', data[[region_var]])
  
  # Define the four regions of interest
  states_of_interest <- c('NSW', 'Qld', 'WA', 'SA', 'Vic')
 
  # Filter data to include only the specified regions
  data <- data %>% filter(State %in% states_of_interest)
  
  # Group data by region
  states <- unique(data[[states_var]])
  numStates <- length(states)
  
  # Initialize result storage
  results <- list()
  FE_estimates <- matrix(list(), nrow = 10, ncol = numRegions) # To store CCE estimates for all regions
  FE_estimates <- matrix(NA, nrow = 10, ncol = numRegions)   # To store CCEMG estimates for each region
  #gdp_percapita_data <- list()                                  # To store GDP per capita data for each region
  
  for (r in 1:numRegions) {
    stateName <- states[r]
    sanitizedstateName <- make.names(stateName, unique = TRUE)
    
    stateData <- data %>% filter(!!sym(states_var) == stateName)
    
    # Get unique years and countries for reshaping
    years <- unique(stateData[[year]])
    statareas <- unique(stateData[[statistical_area]])
    T <- length(years)
    N <- length(statareas)
    
    # Initialize matrices
    Y <- matrix(NA, nrow = T, ncol = N)
    X <- array(NA, dim = c(T, N, 10))
    
    for (i in 1:N) {
      statareaData <- staData %>% filter(!!sym(country_var) == countries[i])
      
      for (t in 1:T) {
        yearData <- countryData %>% filter(!!sym(year_var) == years[t])
        if (nrow(yearData) > 0) {
          if (nrow(yearData) == 1) { # Ensure only a single row is selected
            Y[t, i] <- yearData[[gdp_percapita_var]]
            X[t, i, 1] <- yearData[[labour_force_var]]
            X[t, i, 2] <- yearData[[capital_var]]
            X[t, i, 3] <- yearData[[human_capital_var]]
            X[t, i, 4] <- yearData[[RD_stock_var]]
            X[t, i, 5] <- yearData[[RD_spillovers_var]]
            X[t, i, 6] <- yearData[[infrastructure_var]]
            X[t, i, 7] <- yearData[[urbanization_var]]
          } else {
            warning("Multiple rows found for country: ", countries[i], " and year: ", years[t])
          }
        }
      }
    }
    
    # Compute cross-sectional averages for additional covariates
    cross_sectional_averages <- rowMeans(Y, na.rm = TRUE)
    cross_sectional_averages_lf <- rowMeans(X[, , 1], na.rm = TRUE)
    cross_sectional_averages_hc <- rowMeans(X[, , 3], na.rm = TRUE)
    
    # Add cross-sectional averages to X matrix as new covariates
    for (t in 1:T) {
      X[t, , 8] <- cross_sectional_averages[t]
      X[t, , 9] <- cross_sectional_averages_lf[t]
      X[t, , 10] <- cross_sectional_averages_hc[t]
    }
    
    # Perform CCE and CCEMG estimation
    estimation_result <- CCE_estimator(Y, X, T, N)
    
    # Store results
    results[[sanitizedRegionName]] <- list(
      cce_params = estimation_result$cce_params,
      ccemg_params = estimation_result$ccemg_params,
      ccemg_sd = estimation_result$ccemg_sd / sqrt(N * T),
      ccemg_sdpes = estimation_result$ccemg_sdpes / sqrt(N),
      factors = estimation_result$factors,
      var_ccemg = estimation_result$var_ccemg / (N * T),
      varpes_ccemg = estimation_result$varpes_ccemg / N
    )
    
    # Calculate 95% confidence intervals for CCEMG parameters using both methods
    ci_lower <- estimation_result$ccemg_params - 1.96 * estimation_result$ccemg_sd / sqrt(N * T)
    ci_upper <- estimation_result$ccemg_params + 1.96 * estimation_result$ccemg_sd / sqrt(N * T)
    results[[sanitizedRegionName]]$ccemg_ci <- cbind(ci_lower, ci_upper)
    
    ci_lower_pes <- estimation_result$ccemg_params - 1.96 * estimation_result$ccemg_sdpes / sqrt(N)
    ci_upper_pes <- estimation_result$ccemg_params + 1.96 * estimation_result$ccemg_sdpes / sqrt(N)
    results[[sanitizedRegionName]]$ccemg_ci_pes <- cbind(ci_lower_pes, ci_upper_pes)
    
    # Conduct hypothesis tests for CCEMG parameters
    t_stats <- estimation_result$ccemg_params / (estimation_result$ccemg_sd / sqrt(N * T))
    t_stats_pes <- estimation_result$ccemg_params / (estimation_result$ccemg_sdpes / sqrt(N))
    p_values <- 2 * (1 - pt(abs(t_stats), N - 1)) # Two-tailed test
    p_values_pes <- 2 * (1 - pt(abs(t_stats_pes), N - 1)) # Two-tailed test
    critical_value <- qt(0.975, N - 1) # 95% confidence level
    
    results[[sanitizedRegionName]]$t_stats <- t_stats
    results[[sanitizedRegionName]]$t_stats_pes <- t_stats_pes
    results[[sanitizedRegionName]]$p_values <- p_values
    results[[sanitizedRegionName]]$p_values_pes <- p_values_pes
    results[[sanitizedRegionName]]$critical_value <- critical_value
    
    # Analyze growth patterns
    gdp_growth_rate <- diff(log(Y)) # Log difference to compute growth rate
    mean_growth_rate <- apply(gdp_growth_rate, 1, mean, na.rm = TRUE)
    std_growth_rate <- apply(gdp_growth_rate, 1, sd, na.rm = TRUE)
    
    results[[sanitizedRegionName]]$gdp_growth_rate <- gdp_growth_rate
    results[[sanitizedRegionName]]$mean_growth_rate <- mean_growth_rate
    results[[sanitizedRegionName]]$std_growth_rate <- std_growth_rate
    
    # Compute summary statistics for CCE and CCEMG
    results[[sanitizedRegionName]]$cce_mean <- apply(estimation_result$cce_params, 1, mean, na.rm = TRUE)
    results[[sanitizedRegionName]]$cce_std <- apply(estimation_result$cce_params, 1, sd, na.rm = TRUE)
    
    # Store CCE and CCEMG estimates for all regions
    for (v in 1:10) {
      cce_estimates[v, r] <- list(estimation_result$cce_params[v, ])
      ccemg_estimates[v, r] <- estimation_result$ccemg_params[v]
    }
    
    # Store GDP per capita data for this region
    gdp_percapita_data[[sanitizedRegionName]] <- rowMeans(Y, na.rm = TRUE)
  }
  
  # Determine y-axis limits for each variable across all regions
  y_limits <- matrix(0, nrow = 10, ncol = 2)
  for (v in 1:10) {
    all_estimates <- unlist(cce_estimates[v, ])
    y_limits[v, ] <- range(all_estimates, na.rm = TRUE)
  }
  
  # Plot scatter plots of CCE estimates for each region in one figure for each variable
  for (v in 1:10) {
    par(mfrow = c(ceiling(numRegions / 2), 2))
    for (r in 1:numRegions) {
      regionName <- regions[r]
      scatter_data <- unlist(cce_estimates[v, r])
      plot(1:length(scatter_data), scatter_data, pch = 19, main = paste('Region:', regionName),
           xlab = 'Country', ylab = paste('CCE Estimate of', variable_names[v]), ylim = y_limits[v, ])
      abline(h = ccemg_estimates[v, r], col = 'red', lty = 2, lwd = 2)
      legend('topright', legend = 'CCEMG Estimate', col = 'red', lty = 2, lwd = 2)
    }
    title(main = paste('CCE Estimates for', variable_names[v], 'in All Regions'))
  }
  
  # Plot GDP per capita over time for all regions in the same plot
  plot(NULL, xlim = range(years), ylim = range(unlist(gdp_percapita_data)), xlab = 'Year', ylab = 'GDP per Capita',
       main = 'GDP per Capita over Time for All Regions')
  for (r in 1:numRegions) {
    regionName <- regions[r]
    sanitizedRegionName <- make.names(regionName, unique = TRUE)
    lines(years, gdp_percapita_data[[sanitizedRegionName]], lwd = 2, col = r)
  }
  legend('topright', legend = regions, col = 1:numRegions, lwd = 2)
  
  # Display results summary
  for (r in 1:numRegions) {
    regionName <- regions[r]
    sanitizedRegionName <- make.names(regionName, unique = TRUE)
    
    cat('Results for', regionName, ':\n')
    cat('CCE Mean Estimates:\n')
    print(results[[sanitizedRegionName]]$cce_mean)
    cat('CCE Standard Deviations:\n')
    print(results[[sanitizedRegionName]]$cce_std)
    cat('CCEMG Mean Estimates:\n')
    print(results[[sanitizedRegionName]]$ccemg_params)
    cat('CCEMG Standard Deviations (using Pesaran method):\n')
    print(results[[sanitizedRegionName]]$ccemg_sd)
    cat('CCEMG Standard Deviations (using alternative method):\n')
    print(results[[sanitizedRegionName]]$ccemg_sdpes)
    cat('Variance of CCEMG Estimates (using Pesaran method):\n')
    print(results[[sanitizedRegionName]]$var_ccemg)
    cat('Variance of CCEMG Estimates (using alternative method):\n')
    print(results[[sanitizedRegionName]]$varpes_ccemg)
    cat('95% Confidence Intervals for CCEMG Estimates:\n')
    print(results[[sanitizedRegionName]]$ccemg_ci)
    cat('95% Confidence Intervals for CCEMG Estimates (Pesaran method):\n')
    print(results[[sanitizedRegionName]]$ccemg_ci_pes)
    cat('T-Statistics for CCEMG Estimates (Proposed method):\n')
    print(results[[sanitizedRegionName]]$t_stats)
    cat('T-Statistics for CCEMG Estimates (Pesaran method):\n')
    print(results[[sanitizedRegionName]]$t_stats_pes)
    cat('P-Values for CCEMG Estimates (Proposed method):\n')
    print(results[[sanitizedRegionName]]$p_values)
    cat('P-Values for CCEMG Estimates (Pesaran method):\n')
    print(results[[sanitizedRegionName]]$p_values_pes)
    cat('Critical Value for t-test:', results[[sanitizedRegionName]]$critical_value, '\n')
    cat('Mean GDP Growth Rate:', results[[sanitizedRegionName]]$mean_growth_rate, '\n')
    cat('Standard Deviation of GDP Growth Rate:', results[[sanitizedRegionName]]$std_growth_rate, '\n\n')
  }
  
  # Analyze and compare the results obtained from Pesaran's method and the proposed method
  compare_results(results, regions, variable_names)
  
  # Plot scatter plot of CCE estimates for "Europe and North America and Oceania" region and "ln_rnd_domestic" variable
  plot_scatter_europe_ln_rnd_domestic(results, cce_estimates, ccemg_estimates, regions, 'Europe and North America and Oceania', RD_stock_var)
}

# Function for CCE estimation using plm package
CCE_estimator <- function(Y, X, T, N) {
  # Convert Y and X to a long-format data frame for plm
  data_long <- data.frame(
    year = rep(1:T, N),
    country = rep(1:N, each = T),
    y = as.vector(Y)
  )
  
  for (j in 1:dim(X)[3]) {
    data_long[, paste0("x", j)] <- as.vector(X[, , j])
  }
  
  # Check for zero variance columns
  zero_var_cols <- sapply(data_long, function(col) var(col, na.rm = TRUE) == 0)
  if (any(zero_var_cols)) {
    warning("Zero variance detected in the following columns: ", paste(names(data_long)[zero_var_cols], collapse = ", "))
    data_long <- data_long[, !zero_var_cols]
  }
  
  # Calculate VIF to check for multicollinearity
  formula_str <- paste("y ~", paste(paste0("x", 1:(ncol(data_long) - 2)), collapse = " + "))
  formula <- as.formula(formula_str)
  vif_values <- vif(lm(formula, data = data_long))
  
  # Remove variables with VIF > 10
  high_vif_cols <- names(vif_values[vif_values > 10])
  if (length(high_vif_cols) > 0) {
    warning("High multicollinearity detected. Removing columns: ", paste(high_vif_cols, collapse = ", "))
    data_long <- data_long[, !names(data_long) %in% high_vif_cols]
  }
  
  # Recalculate the formula after removing high VIF columns
  formula_str <- paste("y ~", paste(names(data_long)[!(names(data_long) %in% c("year", "country", "y"))], collapse = " + "))
  formula <- as.formula(formula_str)
  
  # Check if the design matrix is full rank
  design_matrix <- model.matrix(formula, data_long)
  if (qr(design_matrix)$rank < ncol(design_matrix)) {
    stop("The design matrix is rank deficient, indicating perfect multicollinearity.")
  }
  
  # Fit the CCE model using pcce from plm
  model <- plm::pcce(formula, data = data_long, index = c("country", "year"))
  
  # Extract the CCE parameters
  cce_params <- coef(model)
  
  # Compute CCEMG estimates
  ccemg_params <- rowMeans(cce_params, na.rm = TRUE)
  
  # Variance estimation for CCEMG
  covMG1 <- matrix(0, ncol = length(ccemg_params), nrow = length(ccemg_params))
  covMG4 <- matrix(0, ncol = length(ccemg_params), nrow = length(ccemg_params))
  
  for (i in 1:N) {
    if (is.na(cce_params[1, i])) {
      next
    }
    X_i <- X[, i, ]
    Y_i <- Y[, i]
    mat1i <- cce_params[, i] - ccemg_params
    covMG1 <- covMG1 + (mat1i %*% t(mat1i))
    
    e_i <- Y_i - X_i %*% cce_params[, i]
    x_i <- X_i
    
    m <- floor(T^(1/4))
    sumi <- 0
    for (t in 1:m) {
      h_it <- e_i[t] * x_i[t, ]
      sumit <- sumi + h_it %*% t(h_it)
      sumi <- sumit
    }
    
    sigmai0 <- (1/T) * sumi
    suminew <- 0
    for (j in 1:m) {
      weightedi <- 0
      for (t in (j + 1):m) {
        h_it <- e_i[t] * x_i[t, ]
        h_itj <- e_i[t - j] * x_i[t - j, ]
        sumitj <- suminew + (h_it %*% t(h_itj) + h_itj %*% t(h_it))
      }
      weighteditj <- (1 - (j / (m + 1))) * (1 / T) * sumitj
      weightedi <- weightedi + weighteditj
    }
    
    Siepsilon <- sigmai0 + weightedi
    covMG4 <- covMG4 + solve(t(X_i) %*% X_i / T) %*% Siepsilon %*% solve(t(X_i) %*% X_i / T)
  }
  
  estcovMG <- (T / N) * covMG1 + (1 / N) * covMG4
  estcovMGpes <- (1 / N) * covMG1
  var_ccemg <- diag(estcovMG)
  varpes_ccemg <- diag(estcovMGpes)
  ccemg_sd <- sqrt(var_ccemg)
  ccemg_sdpes <- sqrt(varpes_ccemg)
  
  return(list(
    cce_params = cce_params,
    ccemg_params = ccemg_params,
    ccemg_sd = ccemg_sd,
    ccemg_sdpes = ccemg_sdpes,
    factors = model$factors,
    var_ccemg = var_ccemg,
    varpes_ccemg = varpes_ccemg
  ))
}

plot_scatter_europe_ln_rnd_domestic <- function(results, cce_estimates, ccemg_estimates, regions, target_region, variable_name) {
  # Find the index for the target region
  region_idx <- which(regions == target_region)
  if (length(region_idx) == 0) {
    stop(sprintf('Region "%s" not found.', target_region))
  }
  
  # Extract the necessary data
  sanitizedRegionName <- make.names(target_region, unique = TRUE)
  cce_data <- unlist(cce_estimates[4, region_idx]) # 4 corresponds to 'ln_rnd_domestic'
  ccemg_estimate <- ccemg_estimates[4, region_idx]
  ci_proposed <- results[[sanitizedRegionName]]$ccemg_ci[4, ]
  ci_pesaran <- results[[sanitizedRegionName]]$ccemg_ci_pes[4, ]
  
  # Plot the scatter plot
  plot(1:length(cce_data), cce_data, pch = 19, main = sprintf('CCE Estimates for %s in %s', variable_name, target_region),
       xlab = 'Country', ylab = sprintf('CCE Estimate of %s', variable_name))
  abline(h = ccemg_estimate, col = 'red', lty = 2, lwd = 2)
  abline(h = ci_proposed[1], col = 'blue', lty = 3, lwd = 2)
  abline(h = ci_proposed[2], col = 'blue', lty = 3, lwd = 2)
  abline(h = ci_pesaran[1], col = 'green', lty = 3, lwd = 2)
  abline(h = ci_pesaran[2], col = 'green', lty = 3, lwd = 2)
  legend('topright', legend = c('CCEMG Estimate', '95% CI Proposed Lower', '95% CI Proposed Upper', '95% CI Pesaran Lower', '95% CI Pesaran Upper'),
         col = c('red', 'blue', 'blue', 'green', 'green'), lty = c(2, 3, 3, 3, 3), lwd = 2)
}

compare_results <- function(results, regions, variable_names) {
  numRegions <- length(regions)
  
  # Compare standard deviations and variances
  cat('Comparison of Standard Deviations and Variances:\n')
  for (v in 1:10) {
    cat('Variable:', variable_names[v], '\n')
    for (r in 1:numRegions) {
      regionName <- regions[r]
      sanitizedRegionName <- make.names(regionName, unique = TRUE)
      
      ccemg_sd <- results[[sanitizedRegionName]]$ccemg_sd[v]
      ccemg_sdpes <- results[[sanitizedRegionName]]$ccemg_sdpes[v]
      var_ccemg <- results[[sanitizedRegionName]]$var_ccemg[v]
      varpes_ccemg <- results[[sanitizedRegionName]]$varpes_ccemg[v]
      
      cat('Region:', regionName, '\n')
      cat(sprintf('Standard Deviation (Pesaran): %f, Standard Deviation (Proposed): %f\n', ccemg_sd, ccemg_sdpes))
      cat(sprintf('Variance (Pesaran): %f, Variance (Proposed): %f\n', var_ccemg, varpes_ccemg))
    }
    cat('\n')
  }
  
  # Compare confidence intervals
  cat('Comparison of 95% Confidence Intervals:\n')
  for (v in 1:10) {
    cat('Variable:', variable_names[v], '\n')
    for (r in 1:numRegions) {
      regionName <- regions[r]
      sanitizedRegionName <- make.names(regionName, unique = TRUE)
      
      ci <- results[[sanitizedRegionName]]$ccemg_ci[v, ]
      ci_pes <- results[[sanitizedRegionName]]$ccemg_ci_pes[v, ]
      
      cat('Region:', regionName, '\n')
      cat(sprintf('Confidence Interval (Pesaran): [%f, %f]\n', ci[1], ci[2]))
      cat(sprintf('Confidence Interval (Proposed): [%f, %f]\n', ci_pes[1], ci_pes[2]))
    }
    cat('\n')
  }
  
  # Compare t-statistics
  cat('Comparison of T-Statistics:\n')
  for (v in 1:10) {
    cat('Variable:', variable_names[v], '\n')
    for (r in 1:numRegions) {
      regionName <- regions[r]
      sanitizedRegionName <- make.names(regionName, unique = TRUE)
      
      t_stats <- results[[sanitizedRegionName]]$t_stats[v]
      t_stats_pes <- results[[sanitizedRegionName]]$t_stats_pes[v]
      
      cat('Region:', regionName, '\n')
      cat(sprintf('T-Statistic (Proposed): %f, T-Statistic (Pesaran): %f\n', t_stats, t_stats_pes))
    }
    cat('\n')
  }
}
