function plotVolcano(raw_stat, p_val, stat_type, p_threshold, iv_categories, iv_colormap, dvLabels, top_N, fig_size)
% PLOTVOLCANO - Generates a volcano plot with an adaptive vertical-extension 
%               repulsion algorithm and customizable canvas dimensions.
%
% INPUTS:
%   raw_stat      - Vector of the statistical metric (HR, Beta, or T-value)
%   p_val         - Vector of P-Values
%   stat_type     - String: 'HR', 'Beta', or 'T'
%   p_threshold   - Scalar, significance threshold (e.g., 0.05)
%   iv_categories - Cell/categorical array indicating categories for coloring
%   iv_colormap   - Nx3 RGB matrix mapping to the unique categories
%   dvLabels      - Cell array of strings, variable names for labeling
%   top_N         - Number of top significant variables to dynamically label
%   fig_size      - 1x2 Vector: [width, height] in pixels (e.g., [750, 600])

    % ----------------------------------------------------
    % 1. Extract and Process X/Y Data
    % ----------------------------------------------------
    y_data = -log10(p_val);
    
    switch upper(stat_type)
        case 'HR'
            x_data = log2(raw_stat); 
            x_label_str = ['Effect Size: \bf{log_{2}(' stat_type ')}'];
            is_positive_effect = raw_stat > 1; 
            center_line = 0;
            up_legend_text = 'Significant (HR > 1)';
            down_legend_text = 'Significant (HR < 1)';
        case {'BETA', 'T'}
            x_data = raw_stat; 
            x_label_str = ['Effect Size: \bf{' stat_type '}'];
            is_positive_effect = raw_stat > 0; 
            center_line = 0;
            up_legend_text = ['Significant (' stat_type ' > 0)'];
            down_legend_text = ['Significant (' stat_type ' < 0)'];
        otherwise
            error('Invalid stat_type. Choose from ''HR'', ''Beta'', or ''T''.');
    end
    
    point_sizes = y_data * 15 + 10; 

    % ----------------------------------------------------
    % 2. Define Significance Masks and Map Custom Colors
    % ----------------------------------------------------
    is_significant = p_val < p_threshold;
    cats = categorical(iv_categories);
    unique_cats = categories(cats);
    num_cats = length(unique_cats);
    
    cat_colors = iv_colormap; 
    if num_cats > size(cat_colors, 1)
        warning('Number of categories exceeds iv_colormap rows. Appending default colors.');
        cat_colors = [cat_colors; lines(num_cats - size(cat_colors, 1))];
    end
    grey_color = [0.75, 0.75, 0.75]; 

    % ----------------------------------------------------
    % 3. Initialize Figure with Custom Sizes
    % ----------------------------------------------------
    % Fallback to standard size if fig_size is not provided or empty
    if nargin < 9 || isempty(fig_size)
        fig_size = [550, 450]; 
    end
    
    figure('Color', 'w', 'Position', [150, 150, fig_size(1), fig_size(2)]); 
    hold on; grid on;
    ax = gca; ax.LineWidth = 1.2; ax.FontSize = 11;

    % ----------------------------------------------------
    % 4. Plot Non-Significant Points
    % ----------------------------------------------------
    idx_ns_up = (~is_significant) & is_positive_effect;
    scatter(x_data(idx_ns_up), y_data(idx_ns_up), point_sizes(idx_ns_up), ...
        'Marker', '^', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', grey_color, ...
        'MarkerFaceAlpha', 0.4, 'HandleVisibility', 'off');
    
    idx_ns_down = (~is_significant) & (~is_positive_effect);
    scatter(x_data(idx_ns_down), y_data(idx_ns_down), point_sizes(idx_ns_down), ...
        'Marker', 'v', 'MarkerEdgeColor', 'none', 'MarkerFaceColor', grey_color, ...
        'MarkerFaceAlpha', 0.4, 'HandleVisibility', 'off');

    % ----------------------------------------------------
    % 5. Plot Significant Points
    % ----------------------------------------------------
    has_sig_up = false; has_sig_down = false;
    for c = 1:num_cats
        current_cat = unique_cats{c}; current_color = cat_colors(c, :);
        idx_cat_sig = is_significant & (cats == current_cat);
        
        idx_up = idx_cat_sig & is_positive_effect;
        if any(idx_up)
            scatter(x_data(idx_up), y_data(idx_up), point_sizes(idx_up), ...
                'Marker', '^', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', current_color, ...
                'MarkerFaceAlpha', 0.85, 'LineWidth', 0.5, 'HandleVisibility', 'off');
            has_sig_up = true;
        end
        
        idx_down = idx_cat_sig & (~is_positive_effect);
        if any(idx_down)
            scatter(x_data(idx_down), y_data(idx_down), point_sizes(idx_down), ...
                'Marker', 'v', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', current_color, ...
                'MarkerFaceAlpha', 0.85, 'LineWidth', 0.5, 'HandleVisibility', 'off');
            has_sig_down = true;
        end
    end

    % ----------------------------------------------------
    % 6. Add Reference Lines
    % ----------------------------------------------------
    y_limits = ylim;
    plot([center_line center_line], [0 y_limits(2)], '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1, 'HandleVisibility', 'off');
    x_limits = xlim; y_thresh = -log10(p_threshold);
    plot([x_limits(1) x_limits(2)], [y_thresh y_thresh], 'k:', 'LineWidth', 1.2, 'HandleVisibility', 'off');
    text(x_limits(1) + 0.02*(x_limits(2)-x_limits(1)), y_thresh + 0.2, ['p = ' num2str(p_threshold)], 'Color', 'r', 'FontSize', 12, 'FontWeight', 'bold');

   % ----------------------------------------------------
    % 7. Text Labeling with Force-Repulsion Anti-Overlap Algorithm
    % ----------------------------------------------------
    sig_indices = find(is_significant);
    [~, sort_idx] = sort(p_val(sig_indices), 'ascend');
    sorted_sig_indices = sig_indices(sort_idx);
    num_to_label = min(top_N, length(sorted_sig_indices));
    
    if num_to_label > 0
        % Store original point coordinates
        orig_x = x_data(sorted_sig_indices(1:num_to_label));
        orig_y = y_data(sorted_sig_indices(1:num_to_label));
        labels_to_plot = dvLabels(sorted_sig_indices(1:num_to_label));
        
        x_range = x_limits(2) - x_limits(1);
        y_range = y_limits(2) - 0;
        
        % Initialize text positions slightly offset from original points
        text_x = orig_x;
        text_y = orig_y + 0.04 * y_range; % Push up initially
        
        % Map category colors for each specific top label
        label_colors = zeros(num_to_label, 3);
        for i = 1:num_to_label
            target_idx = sorted_sig_indices(i);
            current_cat_val = cats(target_idx);
            cat_idx = find(strcmp(unique_cats, string(current_cat_val)));
            if ~isempty(cat_idx)
                label_colors(i, :) = cat_colors(cat_idx, :);
            else
                label_colors(i, :) = [0 0 0]; % Fallback to black
            end
        end
        
        % Force-directed simulation parameters (Your original framework)
        iterations = 200; 
        repulsion_strength = 0.2 * x_range; 
        min_y_dist = 0.06 * y_range; % Minimum Y distance allowed between texts
        
        % Run repulsion iterations
        for iter = 1:iterations
            for i = 1:num_to_label
                % Dynamically adjust minimum X distance based on text length to prevent overlap
                str_len = length(labels_to_plot{i});
                min_x_dist = (0.04 + 0.008 * str_len) * x_range;
                
                for j = 1:num_to_label
                    if i == j, continue; end
                    
                    % Calculate distance between text bounding approximations
                    dx = text_x(i) - text_x(j);
                    dy = text_y(i) - text_y(j);
                    
                    % If texts are too close, apply repulsive force
                    if abs(dx) < min_x_dist && abs(dy) < min_y_dist
                        % Avoid division by zero
                        if dx == 0, dx = randn * 0.01; end
                        if dy == 0, dy = randn * 0.01; end
                        
                        % Push them apart (Your original linear displacement step)
                        text_x(i) = text_x(i) + sign(dx) * repulsion_strength * 0.1;
                        text_y(i) = text_y(i) + sign(dy) * min_y_dist * 0.1;
                    end
                end
                
                % Boundary check: Prevent labels from flying out or overlapping with the Y-axis
                % Reserve a 6% horizontal safety margin on both left and right sides
                left_bound = x_limits(1) + 0.06 * x_range;
                right_bound = x_limits(2) - 0.06 * x_range;
                
                text_x(i) = max(left_bound, min(right_bound, text_x(i)));
                text_y(i) = max(y_thresh, min(y_limits(2) * 1.05, text_y(i)));
            end
        end
        
        % Draw the optimized labels and connector lines
        for i = 1:num_to_label
            % Draw a colored leader line from original dot to moved text
            plot([orig_x(i), text_x(i)], [orig_y(i), text_y(i)], '-', ...
                'Color', [0.4 0.4 0.4], 'LineWidth', 0.8, 'HandleVisibility', 'off');
            
            % Determine alignment based on its position relative to original dot
            if text_x(i) >= orig_x(i)
                ha = 'left';
            else
                ha = 'right';
            end
            
            % Plot text with matched colors and transparent background box
            text(text_x(i), text_y(i), labels_to_plot{i}, ...
                'FontSize', 10, ...                    
                'FontWeight', 'bold', ...
                'Color', label_colors(i, :), ...      % Matched font color
                'HorizontalAlignment', ha, ...
                'VerticalAlignment', 'middle', ...
                'BackgroundColor', [1 1 1 0.85], ...  % 85% translucent white background box
                'EdgeColor', label_colors(i, :), ...  % Matched border color around text box
                'Margin', 1.5, ...                      
                'Interpreter', 'none');
        end
    end



    % ----------------------------------------------------
    % 8. Legend and Labels
    % ----------------------------------------------------
    legend_handles = []; legend_labels = {};
    if has_sig_up
        h_dummy_up = scatter(NaN, NaN, 60, 'Marker', '^', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0.3 0.3 0.3], 'LineWidth', 0.8);
        legend_handles = [legend_handles, h_dummy_up]; legend_labels = [legend_labels, {up_legend_text}];
    end
    if has_sig_down
        h_dummy_down = scatter(NaN, NaN, 60, 'Marker', 'v', 'MarkerEdgeColor', 'k', 'MarkerFaceColor', [0.3 0.3 0.3], 'LineWidth', 0.8);
        legend_handles = [legend_handles, h_dummy_down]; legend_labels = [legend_labels, {down_legend_text}];
    end
    if ~isempty(legend_handles), legend(legend_handles, legend_labels, 'Location', 'best'); end

    xlabel(x_label_str, 'FontSize', 12);
    ylabel('Significance: \bf{-log_{10}(P-Value)}', 'FontSize', 12);
    title(['Volcano Plot of ' stat_type ' Results'], 'FontSize', 14, 'FontWeight', 'bold');
    
    % Dynamically expand Y-limits if labels extend above the standard viewport
    if num_to_label > 0
        ylim([0, max(y_limits(2), max(text_y)*1.15)]);
    end
    hold off;
end