function [ c_k_plus_1 ] = proximal_mapping( S_upper, proximal_mapping_term )
%PROXIMAL_MAPPING Summary of this function goes here
%   Detailed explanation goes here

lower_bounds = zeros(length(proximal_mapping_term),1);
% lower_bounds = ones(length(proximal_mapping_term),1)*100;

upper_bounds = zeros(length(proximal_mapping_term),1) + S_upper;

max_operation = [proximal_mapping_term lower_bounds];

proximal_mapping_term = max(max_operation,[],2);

min_operation = [proximal_mapping_term upper_bounds];

proximal_mapping_term = min(min_operation,[],2);

projection_counter = 0;

% S_upper

while sum(proximal_mapping_term) - S_upper > 1e-5
    
    positive_indices = logical(proximal_mapping_term);
    
    lambda = (sum(proximal_mapping_term)-S_upper)/(length(proximal_mapping_term) - length(find(~positive_indices)));
    
    proximal_mapping_term(positive_indices) = proximal_mapping_term(positive_indices) - lambda;
    
    max_operation = [proximal_mapping_term lower_bounds];
    
    proximal_mapping_term = max(max_operation,[],2);
    
    min_operation = [proximal_mapping_term upper_bounds];
    
    proximal_mapping_term = min(min_operation,[],2);
    
    projection_counter = projection_counter + 1;
    
end

c_k_plus_1 = proximal_mapping_term;

end

