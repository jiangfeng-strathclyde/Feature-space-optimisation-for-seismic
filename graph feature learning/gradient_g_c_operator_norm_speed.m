function [ grad_g_c, g_c, L_constant ] = gradient_g_c_operator_norm_speed( c_k, feature, FD, x )
%G_C Summary of this function goes here
%   Detailed explanation goes here

grad_g_c = zeros(size(feature,2),1);
grad_grad_g_c = zeros(size(feature,2),size(feature,2));

g_c = 0;
x = sign(x);

feature_weight = zeros(size(feature,1));
feature_weight_lower = zeros(size(feature,1));
x_r_c = zeros(size(feature,1));


for i_ggc = 1:size(feature,2)
    
    if i_ggc == 1
        
        x_r_c = (x - x.').^2;

        
        
        for rc_i = 1:size(feature,2)
            
            feature_weight = feature_weight - c_k(rc_i)*FD{rc_i};
            
            feature_weight_lower = feature_weight_lower - 1 *FD{rc_i};
            
        end
     
        
        
        W = exp(feature_weight);
        
        W(W == diag(W)) = 0;

        W_lower = exp(feature_weight_lower);
        
        W_lower(W_lower == diag(W_lower)) = 0;
        
    end
    
    grad_g_c_temp = W .* FD{i_ggc} .* x_r_c;
    
    grad_g_c(i_ggc) = - sum(sum(grad_g_c_temp));
    
    grad_g_c_lower_temp = W_lower .* FD{i_ggc} .* x_r_c;
    
    

    for j_ggc = 1:size(feature,2)
        
        grad_grad_g_c_temp = grad_g_c_lower_temp .* FD{j_ggc};
        
        grad_grad_g_c(i_ggc,j_ggc) = sum(sum(grad_grad_g_c_temp));
        
    end
 
    
end

%% Lipschitz constant
L_constant = sqrt(max(eig(grad_grad_g_c' * grad_grad_g_c))); 

end

















