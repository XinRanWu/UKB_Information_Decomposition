function [result_table] = mean_res_syn_intra_extra_globe(info_table,id_net,redundancy,synergy)
    n_subject = size(redundancy,3);
    n_region = size(redundancy,1);
    id_net_right = id_net(id_net <= 180);
    id_net_left = id_net(id_net > 180);
    id_right_left = {id_net_right;id_net_left};
    
    intra_right_idx = zeros(n_region); 
    intra_right_idx(id_right_left{1},id_right_left{1}) = 1;

    intra_left_idx = zeros(n_region); 
    intra_left_idx(id_right_left{2},id_right_left{2}) = 1;

    intra_left_right_idx = zeros(n_region);
    intra_left_right_idx(id_right_left{1},id_right_left{1}) = 1;
    intra_left_right_idx(id_right_left{2},id_right_left{2}) = 1;

    extra_left_right_idx = zeros(n_region);
    extra_left_right_idx(id_right_left{1},id_right_left{2}) = 1;

    intra_right_red = zeros(n_subject,1);
    intra_right_syn = zeros(n_subject,1);
    intra_left_red = zeros(n_subject,1);
    intra_left_syn = zeros(n_subject,1);
    intra_left_right_red = zeros(n_subject,1);
    intra_left_right_syn = zeros(n_subject,1);
    extra_left_right_red = zeros(n_subject,1);
    extra_left_right_syn = zeros(n_subject,1);
    for i = 1:n_subject
        red_mat = redundancy(:,:,i);
        syn_mat = synergy(:,:,i);
    
        intra_right_red_mat = red_mat.*intra_right_idx;
        intra_right_syn_mat = syn_mat.*intra_right_idx;
    
        intra_left_red_mat = red_mat.*intra_left_idx;
        intra_left_syn_mat = syn_mat.*intra_left_idx;
    
        intra_left_right_red_mat = red_mat.*intra_left_right_idx;
        intra_left_right_syn_mat = syn_mat.*intra_left_right_idx;
    
        extra_left_right_red_mat = red_mat.*extra_left_right_idx;
        extra_left_right_syn_mat = syn_mat.*extra_left_right_idx;

        intra_right_red(i) = mean(intra_right_red_mat(intra_right_red_mat~=0),'all');
        intra_right_syn(i) = mean(intra_right_syn_mat(intra_right_syn_mat~=0),'all');
        intra_left_red(i) = mean(intra_left_red_mat(intra_left_red_mat~=0),'all');
        intra_left_syn(i) = mean(intra_left_syn_mat(intra_left_syn_mat~=0),'all');
        intra_left_right_red(i) = mean(intra_left_right_red_mat(intra_left_right_red_mat~=0),'all');
        intra_left_right_syn(i) = mean(intra_left_right_syn_mat(intra_left_right_syn_mat~=0),'all');
        extra_left_right_red(i) = mean(extra_left_right_red_mat(extra_left_right_red_mat~=0),'all');
        extra_left_right_syn(i) = mean(extra_left_right_syn_mat(extra_left_right_syn_mat~=0),'all');
    end

    result_table = addvars(info_table, intra_right_red,intra_right_syn,...,
        intra_left_red,intra_left_syn,intra_left_right_red,intra_left_right_syn,extra_left_right_red,...,
        extra_left_right_syn,'NewVariableNames', {'intra_right_red','intra_right_syn',...,
        'intra_left_red','intra_left_syn','intra_left_right_red','intra_left_right_syn',...,
        'extra_left_right_red','extra_left_right_syn'});
end

