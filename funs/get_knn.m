function [B] = get_knn(X,centers,k)
D = L2_distance_1(X', centers');
[~, idx] = sort(D, 2); % sort each row
n=size(X,1);
m=size(centers,1);

B = zeros(n,m);
for ii = 1:n
    id = idx(ii,1:k+1);
    di = D(ii, id);
    B(ii,id) = (di(k+1)-di)/(k*di(k+1)-sum(di(1:k))+eps);
end
end