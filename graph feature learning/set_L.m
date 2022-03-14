function [ L, W ] = set_L( feature, FD, c_k )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

W = zeros(size(feature,1));

for j = 1:size(feature,2)
    
    W = W + c_k(j) * FD{j};
    
end

W = exp(-W);


% hist(W,10000)


% W(W == diag(W)) = 0;

%  W = W - min(W(:));
%  W = W/max(W(:));
 W(W == diag(W)) = 0; 


L = diag(sum(W))-W;

end

