function [A] = Threedcorr(X,Y)
% This function calculates the correlation values between two 3D arrays of the same dimensions. 
% Arrays X and Y have dimensions [A,A,B]. 
% Firstly, the upper triangular part of the first two dimensions is extracted and flattened into a one-dimensional vector. 
% Considering the third dimension, the intermediate result has dimensions [A(A - 1)/2, B]. 
% Then, Pearson correlation coefficients are computed row-wise for these two intermediate matrices. 
% Finally, a vector of dimensions [B,1] is obtained.

X = Triu_three(X);
Y = Triu_three(Y);
reshaped_X = reshape(X, [], size(X, 3));
reshaped_triu_X = reshaped_X(~all(reshaped_X == 0, 2), :);
reshaped_Y = reshape(Y, [], size(Y, 3));
reshaped_triu_Y = reshaped_Y(~all(reshaped_X == 0, 2), :);

% Here,The columns of reshaped_triu_X correspond to the third dimension of 
% the original X, and the rows correspond to the upper triangular elements 
% of the first two dimensions.

A = diag(corr([reshaped_triu_X,reshaped_triu_Y]), size(reshaped_triu_X,2) );


% ************************************************
% Obtain the upper triangular part of the first two dimensions of a three-dimensional matrix.
% ************************************************
function [M] = Triu_three(Z)
    for i = 1:size(Z,3)
        Z(:,:,i) = triu(Z(:,:,i), 1);
    end
M = Z;
end

end