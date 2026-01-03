function array = network2array(conn, have_diagonal)
%NETWORK2ARRAY 提取矩阵的上三角
%   Detailed explanation goes here

n = length(conn);
if have_diagonal
    array = zeros([1 n*(n+1)/2]);
    k=1;
    for i=1:n
        for j=i:n
            array(k) = conn(i,j);
            k=k+1;
        end
    end
else
    array = zeros([1 n*(n-1)/2]);
    k=1;
    for i=1:n
        for j=i+1:n
            array(k) = conn(i,j);
            k=k+1;
        end
    end
end

end

