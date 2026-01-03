function [A] = Matrix_reduce(X,id_cell)
% This function is designed to compute reduction of a two-dimensional matrix 
% given index groups; X is the two-dimensional symmetric matrix to be operated on, 
% id_cell is the cell containing grouping IDs.

N_group = size(id_cell,1);
A = zeros(N_group); 
for i = 1: N_group
    for j = 1: N_group
        A(i,j) = mean(X(id_cell{i},id_cell{j}),'all');
    end
end

end

