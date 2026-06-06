function [dia] = icatb_mat2dia3d(mat)
dia = zeros(size(mat,3),size(mat(:,:,1),1));
id = find(eye(size(mat(:,:,1))));
for i = 1:size(mat,3)
    m = mat(:,:,i);dia(i,:) = m(id);
    % disp(num2str(i));
end
end