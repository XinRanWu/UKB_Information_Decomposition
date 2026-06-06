function plot_interEff(data, yVar, xVar, mVar, covariates, lineColor)
% plot_interEff Plots the interaction effect of two continuous variables with accurate CIs.
% Supports both fitlm and fitlme (automatically handles Centre_0_0 as a random effect).

    %% Step 0: Check if all variables exist in the table
    allVars = [{yVar, xVar, mVar}, covariates];
    missingVars = allVars(~ismember(allVars, data.Properties.VariableNames));
    if ~isempty(missingVars)
        error('The following variables were not found in your table: %s', strjoin(missingVars, ', '));
    end

    %% Step 1: Regress out covariates if they exist (get residuals)
    if ~isempty(covariates)
        % Check if Centre_0_0 is in covariates to decide between fitlm and fitlme
        hasCentre = any(strcmp('Centre_0_0', covariates));
        
        if hasCentre
            % Remove Centre_0_0 from fixed effects and add it as a random effect
            fixedCovs = covariates(~strcmp('Centre_0_0', covariates));
            covStr = strjoin(fixedCovs, ' + ');
            residualFormula = [yVar, ' ~ ', covStr, ' + (1|Centre_0_0)'];
            fprintf('Regressing out covariates using MIXED MODEL (fitlme).\nFormula: %s\n', residualFormula);
            
            residualModel = fitlme(data, residualFormula);
        else
            % Standard linear model if no random effect is needed
            covStr = strjoin(covariates, ' + ');
            residualFormula = [yVar, ' ~ ', covStr];
            fprintf('Regressing out covariates using LINEAR MODEL (fitlm).\nFormula: %s\n', residualFormula);
            
            residualModel = fitlm(data, residualFormula);
        end
        
        % Replace original Y with its residuals
        data.Y_analyzed = residuals(residualModel);
    else
        data.Y_analyzed = data.(yVar);
    end

    %% Step 2: Construct and fit the interaction model
    % Rename columns uniformly to prevent formula parsing issues
    data.X_analyzed = data.(xVar);
    data.M_analyzed = data.(mVar);
    
    interactionFormula = 'Y_analyzed ~ X_analyzed * M_analyzed';
    model = fitlm(data, interactionFormula);
    
    % Extract coefficients and covariance matrix
    beta = model.Coefficients.Estimate;
    covMat = model.CoefficientCovariance;

    %% Step 3: Define Moderator (M) levels (Mean ¡À 1 SD)
    M_val = data.M_analyzed;
    M_mean = mean(M_val, 'omitnan');
    M_std  = std(M_val, 'omitnan');
    
    M_levels = [M_mean - M_std, M_mean, M_mean + M_std];
    M_labels = {'Low (-1 SD)', 'Mean', 'High (+1 SD)'};
    lineStyles = {'-.', '--', '-'};

    %% Step 4: Prepare X-axis range for predictions
    X_val = data.X_analyzed;
    x_range = linspace(min(X_val), max(X_val), 100)';
    
    % Initialize plot properties
    hold on;
    grid on;

    %% Step 5: Calculate and plot curves with CIs for each level
    for i = 1:3
        m_curr = M_levels(i);
        
        % Construct design matrix: [1, X, M, X*M]
        X_design = [ones(size(x_range)), x_range, ones(size(x_range))*m_curr, x_range*m_curr];
        
        % Calculate fitted Y values
        y_fit = X_design * beta;
        
        % Exact calculation of standard errors
        y_se = sqrt(sum((X_design * covMat) .* X_design, 2));
        
        % Calculate 95% Confidence Intervals
        ci_low = y_fit - 1.96 * y_se;
        ci_high = y_fit + 1.96 * y_se;
        
        % 1. Plot the confidence interval shading
        fill([x_range; flipud(x_range)], [ci_low; flipud(ci_high)], ...
             lineColor, 'FaceAlpha', 0.08, 'EdgeColor', 'none', 'HandleVisibility', 'off');
         
        % 2. Plot the main regression line
        plot(x_range, y_fit, lineStyles{i}, 'Color', lineColor, 'LineWidth', 2.5, ...
             'DisplayName', sprintf('%s = %.2f (%s)', mVar, m_curr, M_labels{i}));
    end

    %% Step 6: Formatting and aesthetics
    xlabel(xVar, 'Interpreter', 'none');
    ylabel(['Residuals of ', yVar], 'Interpreter', 'none');
    title(sprintf('Interaction Effect on %s', yVar), 'Interpreter', 'none');
    legend('Location', 'best');
    hold off;
end