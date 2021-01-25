% Fornino Manera (2019)
% 
% DATE: January 24, 2021
%
% Project: Automation and the Future of Work: Assessing the Role of Labor
%          Flexibility

function ndx = sub2ind_vec(siz, vectors)
% INVERTS IND2SUB
% INPUT:
% vectors: a matrix of vectors containing the indices of a grid of points like
% [s_3 s_2 s_1]
% siz: a vector giving the size [n_3 n_2 n_1] of the grid of points
% n_1 * n_2 * n_3, where n_1 is the most external state so that e.g.
% position 2 will correspond to s_1(1), s_2(1), s_3(2)
% position n_2 + 1 will be: s_1(1), s_2(2), s_3(1)
% position n_2 * n_3 + 1 : s_1(2), s_2(1), s_3(1)
% 
% OUTPUT:
% ndx: the indices on the grid n_1 * n_2 * n_3 corresponding to the above


dim = length(siz);
ndx = vectors(:,1);
k = cumprod(siz);
for i = 2:dim
    ndx = ndx + k(i-1) * (vectors(:,i) - 1);
end

end
