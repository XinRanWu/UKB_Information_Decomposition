function [net_mean_mat] = NetMeanMatrix3D(mat, net_idx)
    net_idx = net_idx(:); 
    net_num = unique(net_idx);
    num_nets = length(net_num);
    num_slices = size(mat, 3);
    net_mean_mat = zeros(num_nets, num_nets, num_slices);
    net_masks = cell(num_nets, 1);
    for i = 1:num_nets
        net_masks{i} = (net_idx == net_num(i));
    end
    for i = 1:num_nets
        mask_i = net_masks{i};
        for j = 1:num_nets
            mask_j = net_masks{j};
            if i == j
                n_nodes = sum(mask_i);
                if n_nodes <= 1
                    net_mean_mat(i,j,:) = 0;
                    continue;
                end
 
                sub_mat = mat(mask_i, mask_i, :);
                lower_tri_idx = find(tril(ones(n_nodes), -1));
                sub_mat_flat = reshape(sub_mat, n_nodes * n_nodes, num_slices);
                net_mean_mat(i,j,:) = mean(sub_mat_flat(lower_tri_idx, :), 1);
                
            else
                mean_dim1 = mean(mat(mask_i, mask_j, :), 1);
                net_mean_mat(i,j,:) = mean(mean_dim1, 2);    
            end
        end
    end
end