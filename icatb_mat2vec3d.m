function [vec] = icatb_mat2vec3d(mat)
vec = zeros(size(mat,3),length(icatb_mat2vec(mat(:,:,1))));
for i = 1:size(mat,3)
    vec(i,:) = icatb_mat2vec(mat(:,:,i));
    % disp(num2str(i));
end
end