% Fornino Manera (2019)
% 
% DATE: January 24, 2021
%
% Project: Automation and the Future of Work: Assessing the Role of Labor
%          Flexibility

function vector = ind2sub_vec(siz, ndx)
% INPUT:
% siz: a vector containing the length of a grid of points like
% [n_3 n_2 n_1]
% n_1 * n_2 * n_3, where n_1 is the most external state so that e.g.
% position 2 will correspond to s_1(1), s_2(1), s_3(2)
% position n_2 + 1 will be: s_1(1), s_2(2), s_3(1)
% position n_2 * n_3 + 1 : s_1(2), s_2(1), s_3(1)
% ndx: the indices to query for
%
% OUTPUT:
% vector: a (length(ndx))-(length(siz)) matrix of indices for s_3 s_2 s_1 


dim = length(siz);
queries = length(ndx);
vector = zeros(queries,dim);
k = cumprod(siz);
for i = dim:-1:2
    remainder = rem(ndx - 1, k(i-1)) + 1;
    vector(:,i) = (ndx - remainder) / k(i-1) + 1;
    ndx = remainder;
end
vector(:,1) = ndx;

end
