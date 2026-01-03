function [mean_array] = mean_extra_net(mat,id_net)
    n_subject = size(mat,3);
    n_region = size(mat,1);
    extra_net_idx = zeros(n_region);
    extra_net_idx(id_net,:) = 1;
    extra_net_idx(id_net,id_net) = 0;
    mean_array = zeros(n_subject,1);
    for i = 1:n_subject
        data = mat(:,:,i);
        extra_net_mat = data.*extra_net_idx;
        mean_array(i) = mean(extra_net_mat(extra_net_mat ~= 0),'all');
    end
end

