function [results_table] = CoxRegression(input_table, indep_vars, dep_vars, covariate_vars)
% BATCHCOXSURVIVAL - Batch Cox Proportional Hazards Regression using coxphfit
%   INPUTS:
%       input_table        - Data table (e.g., UKB_BrainXY)
%       indep_vars         - Cell array of independent variable names
%       dep_vars           - Cell array of survival/event time columns
%       covariate_vars     - Cell array of covariate names (Sex, Age, etc.)
%       max_follow_up_time - Scalar or vector representing the maximum follow-up 
%                            time for censored individuals (NaN cases)
%   OUTPUT:
%       results_table      - Table containing batch Cox regression results

warning('off');
results = {};
row = 1;

for d = 1:length(dep_vars)
    dep_col = dep_vars{d};
    max_follow_up_time = max(dep_col);
    for i = 1:length(indep_vars)
        iv_col = indep_vars{i};
        fprintf('Analyzing: %s | %s\n', dep_col, iv_col);
        
        time_raw    = input_table.(dep_col);
        iv_data     = input_table.(iv_col);
        
        % Extract covariate data
        cov_data    = [];
        for c = 1:length(covariate_vars)
            cov_data = [cov_data, input_table.(covariate_vars{c})];
        end
        
        % ----------------------------------------------------
        % 1. Data Cleaning: Keep samples with complete IV and Covariates
        % ----------------------------------------------------
        valid_mask = ~isnan(iv_data);
        for c = 1:size(cov_data,2)
            valid_mask = valid_mask & ~isnan(cov_data(:,c));
        end
        
        time_filtered = time_raw(valid_mask);
        iv_filtered   = iv_data(valid_mask);
        cov_filtered  = cov_data(valid_mask,:);
        
        % ----------------------------------------------------
        % 2. Define Survival Time and Censoring Status
        %    MATLAB convention: 1 = Censored, 0 = Event occurred
        % ----------------------------------------------------
        final_time = zeros(size(time_filtered));
        censoring  = zeros(size(time_filtered)); 
        
        for k = 1:length(time_filtered)
            if isnan(time_filtered(k))
                % NaN means no event occurred during follow-up (Censored)
                censoring(k) = 1; 
                % Assign maximum follow-up time
                if isscalar(max_follow_up_time)
                    final_time(k) = max_follow_up_time;
                else
                    max_follow_up_filtered = max_follow_up_time(valid_mask);
                    final_time(k) = max_follow_up_filtered(k);
                end
            else
                % Non-NaN means event occurred
                censoring(k) = 0;
                final_time(k) = abs(time_filtered(k));
            end
        end
        
        % ----------------------------------------------------
        % 3. Exclude Invalid Survival Time (Time <= 0)
        % ----------------------------------------------------
        keep = final_time > 1e-6;
        if sum(keep) < 5 || sum(censoring(keep) == 0) < 3
            % Skip if sample size or event count is too small
            fprintf('  -> Sample or event size too small, skipped.\n');
            continue;
        end
        
        X_final = [iv_filtered(keep), cov_filtered(keep,:)];
        T_final = final_time(keep);
        C_final = censoring(keep);
        
        % ----------------------------------------------------
        % 4. Model Fitting with Error Handling
        % ----------------------------------------------------
        try
            [b, logL, ~, stats] = coxphfit(X_final, T_final, 'Censoring', C_final);
            
            beta    = b(1);
            se      = stats.se(1);
            p       = stats.p(1);
            hr      = exp(beta);
            hr_low  = exp(beta - 1.96*se);
            hr_high = exp(beta + 1.96*se);
            n       = sum(keep);
            n_event = sum(C_final == 0); 
            
            % Store statistics into results cell array
            results{row, 1}  = dep_col;
            results{row, 2}  = iv_col;
            results{row, 3}  = beta;
            results{row, 4}  = se;
            results{row, 5}  = hr;
            results{row, 6}  = hr_low;
            results{row, 7}  = hr_high;
            results{row, 8}  = p;
            results{row, 9}  = logL;
            results{row, 10} = n;
            results{row, 11} = n_event;
            row = row + 1;
        catch ME
            % Prevent batch crash if optimization fails to converge
            fprintf('  -> Fit failed for %s: %s\n', iv_col, ME.message);
        end
    end
end

% ----------------------------------------------------
% 5. Build Result Table
% ----------------------------------------------------
if isempty(results)
    results_table = table();
    fprintf('\nNo valid results generated.\n');
else
    results_table = cell2table(results, ...
        'VariableNames', {'Survival_Var','Independent_Var',...
        'Coefficient','SE','HR','HR_95CI_Lower','HR_95CI_Upper',...
        'P_Value','Log_Likelihood','Sample_Size','Event_Size'});
    fprintf('\nAnalysis done!\n');
end
end