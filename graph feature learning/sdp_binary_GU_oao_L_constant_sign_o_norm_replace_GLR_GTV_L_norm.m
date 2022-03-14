function [ x_valid, class_SDP_temp, SDP_error,ck, whole_final_ck] = ...
    sdp_binary_GU_oao_L_constant_sign_o_norm_replace_GLR_GTV_L_norm( class_test, ...
    c_k_initially_set, tol_set, tol_set_pg, S_lower, S_upper, ...
    feature_train_test, initial_label, initial_label_index, class_train_test, alpha ,training_label,sigma)
global labelednum unlabelednum
%SDP_BINARY_GU_oao Summary of this function goes here
%   Detailed explanation goes here
% classname={'class1.1','class1.2','class2.2','class2.3','class2.4','class3.1','class3.2','class3.4','class9'};
%% set L (graph Laplacian (based on features that are drived from body part trajectories))
tol = 1e+04;
counter = 0;
GT_obj_all = 0;
SDP_obj_all = 0;
SDP_error = 0;
% sigma=1;
trade_off=1;
whole_final_ck=zeros(119,1000);
ck_iteration=0;

FD = cell(size(feature_train_test,2),1);
diff = cell(size(feature_train_test,2),1);
for j = 1:size(feature_train_test,2)
  
    FD{j} = (feature_train_test(:,j) - feature_train_test(:,j).').^2;
    
end
c_k_initially_set=1/sigma^2;
while tol>tol_set
    
    if counter == 0
        
        %% ==== SET INITIAL FEATURE WEIGHTS ALL TO BE 1 =====
        c_k = zeros(size(feature_train_test,2),1)+c_k_initially_set;

        %% =================================================      
        c_k_no_GU = c_k;

        %=========SDP classifier learning==========
        % set graph Laplacian
        [L,W] = set_L( feature_train_test, FD, c_k );
        D=diag(sum(W));
        D_diag=diag(D);
        [z_r,z_c]=find(D_diag==0);
        D(z_r,z_r)=mean(D_diag);
        D_diag=diag(D);    
        L(isnan(L))=0;
        L(isinf(L))=0;
        
        half_D=zeros(size(D));
        half_D(logical(eye(size(half_D))))=D_diag.^(-1/2);
        L_norm=half_D*L*half_D;
        L_norm(isnan(L_norm))=0;
        L_norm(isinf(L_norm))=0;
        
        n = size(class_train_test,1);
        W_normalized = W./max(W(:));
        x_known=[training_label;ones(length(class_train_test)-length(training_label),1)*0];
        [nRows,nCols] = size(W_normalized);
        C_matrix=zeros(nRows,nCols);
        subsetIdx = 1:length(training_label);
        diagonalIdx = (subsetIdx-1) * (nRows + 1) + 1;
        C_matrix(diagonalIdx) = 1;

%% CVX begin to upload label
        
    L=L; % l-matrix
    cvx_begin sdp
        variable nu(n)
        maximize ( -sum(nu) )
        L+ diag(nu) >= 0;
    cvx_end

        x=x_known;
        num_train=length(find(x_known~=0));
        num_test=length(find(x_known==0));
        new_L=L+diag(nu);
        W3=new_L(num_train+1:end,num_train+1:end);
        W1=new_L(1:num_train,1:num_train);
        W2=new_L(1:num_train,num_train+1:end);
        x1=x(1:num_train);
        x2=x(num_train+1:end);

        B=W2'*x1;
        C=x1'*W2;
        new_x=-inv(W3)*(B+C')/2;
        x=[training_label;new_x];
    
    %% cvx end
        
        check_obj_with_GT_label = trace(L*(class_train_test*class_train_test'));
        
        GT_obj_all(counter+1)=check_obj_with_GT_label;
        
        SDP_obj_all(counter+1)=cvx_optval;
        
        GT_no_GU_obj = check_obj_with_GT_label;
        SDP_no_GU_obj = cvx_optval;
        
        class_SDP_no_GU_temp_binary = sign(x);

        class_SDP_no_GU_temp_binary(initial_label_index) = [];
        
        class_SDP_no_GU_temp = zeros(size(class_SDP_no_GU_temp_binary,1),size(class_SDP_no_GU_temp_binary,2));
        class_SDP_no_GU_temp(class_SDP_no_GU_temp_binary==1) = 1;
        class_SDP_no_GU_temp(class_SDP_no_GU_temp_binary==-1) = -1;
        
        x_valid = x;
        %==========================================
        
       diff_label = sign(x_valid(964:end))-class_train_test(964:end);
       conmat=confusionmat(class_train_test(964:end),sign(x_valid(964:end)));

        error_SDP = size(find(diff_label~=0),1)*size(find(diff_label~=0),2)/size(class_test,1);
        
        SDP_error(counter+1)=error_SDP;
        
        SDP_no_GU_error = error_SDP;
        
        disp(['initial classifier accuracy: ']);
        disp(conmat);

        disp(['initial objective: ' num2str(x_valid'*L*x_valid)]);
        
        pause(0.0001);
    end
    
    tol_pg = 1e+4;
    [grad_g_c, ~, Lipschitz_constant] = gradient_g_c_operator_norm_speed( c_k, feature_train_test, FD, x );
    
   step_size =  2/Lipschitz_constant;
    precision = step_size;
    
    counter_pg = 0;
     while tol_pg > tol_set_pg
        
        %% ================GRAPH UPDATE STARTS FROM HERE==============
        
        if counter_pg > 0
            [grad_g_c, ~, ~] = gradient_g_c_operator_norm_speed( c_k, feature_train_test, FD, x );
        end
        %% proximal gradient
        ck_iteration = ck_iteration+1;
        whole_final_ck(:,ck_iteration) = c_k;
        c_k_REAL = c_k;
        
        proximal_mapping_term = c_k_REAL - step_size * grad_g_c;

        c_k_plus_1 = proximal_mapping( S_upper, proximal_mapping_term );
        obj_previous = sign(x)' * L * sign(x);
        [L,W] = set_L( feature_train_test, FD, c_k_plus_1 );
       
        D=diag(sum(W));
        D_diag=diag(D);
        [z_r,z_c]=find(D_diag==0);
        D(z_r,z_r)=mean(D_diag);
        D_diag=diag(D);    
        L(isnan(L))=0;
        L(isinf(L))=0;
        
        half_D=zeros(size(D));
        half_D(logical(eye(size(half_D))))=D_diag.^(-1/2);
        L_norm=half_D*L*half_D;
        L_norm(isnan(L_norm))=0;
        L_norm(isinf(L_norm))=0;
        W_normalized = W./max(W(:));
        obj_current = sign(x)' * L * sign(x);
        tol_pg = abs(obj_current - obj_previous); % update tol_pg

        c_k = c_k_plus_1;
        
        counter_pg = counter_pg + 1;
        disp(tol_pg)
        pause(0.00001);
        
     end
    
    counter = counter + 1;

    %=========SDP classifier learning==========
    % set graph Laplacian
    [L,W] = set_L( feature_train_test, FD, c_k );
    n = size(class_train_test,1);
        D=diag(sum(W));
        D_diag=diag(D);
        [z_r,z_c]=find(D_diag==0);
        D(z_r,z_r)=mean(D_diag);
        D_diag=diag(D);    
        L(isnan(L))=0;
        L(isinf(L))=0;
        
        half_D=zeros(size(D));
        half_D(logical(eye(size(half_D))))=D_diag.^(-1/2);
        L_norm=half_D*L*half_D;
        L_norm(isnan(L_norm))=0;
        L_norm(isinf(L_norm))=0;
        
    n = size(class_train_test,1);
        W_normalized = W./max(W(:));
        x_known=[training_label;ones(length(class_train_test)-length(training_label),1)*0];
        [nRows,nCols] = size(W_normalized);
        C_matrix=zeros(nRows,nCols);
        subsetIdx = 1:length(training_label);
        diagonalIdx = (subsetIdx-1) * (nRows + 1) + 1;
        C_matrix(diagonalIdx) = 1;
        T=isnan(L_norm);
        [row_L_norm,~]=size(L_norm);
    
        if sum(sum(T))~=0
            L_norm=ones(row_L_norm);
        end
        
%% CVX begin
    
L=L; % l-matrix
    cvx_begin sdp
        variable nu(n)
        maximize ( -sum(nu) )
        L+ diag(nu) >= 0;
    cvx_end

        x=x_known;
        num_train=length(find(x_known~=0));
        num_test=length(find(x_known==0));
        new_L=L+diag(nu);
        W3=new_L(num_train+1:end,num_train+1:end);
        W1=new_L(1:num_train,1:num_train);
        W2=new_L(1:num_train,num_train+1:end);
        x1=x(1:num_train);
        x2=x(num_train+1:end);

        B=W2'*x1;
        C=x1'*W2;
        new_x=-inv(W3)*(B+C')/2;
        x=[training_label;new_x];
        
        %% CVX end

    check_obj_with_GT_label = trace(L*(class_train_test*class_train_test'));
    
    GT_obj_all(counter+1)=check_obj_with_GT_label;
    
    SDP_obj_all(counter+1)=cvx_optval;
    
    if isempty(find(isnan(SDP_obj_all(counter+1))==1)) == 1
        x_valid = sign(x);
    else
    end
    x_valid = x;

diff_label = sign(x_valid(labelednum+1:end))-class_train_test(labelednum+1:end);
       conmat=confusionmat(class_train_test(labelednum+1:end),sign(x_valid(labelednum+1:end)));

        error_SDP = size(find(diff_label~=0),1)*size(find(diff_label~=0),2)/size(class_test,1);
        
        SDP_error(counter+1)=error_SDP;
        
        SDP_no_GU_error = error_SDP;
        
        disp(['initial classifier accuracy: ']);
        disp(conmat);

        disp(['initial objective: ' num2str(x_valid'*L*x_valid)]);
        
        pause(0.0001);
    
    %==========================================
    disp(['classifier ' num2str(counter) ' accuracy: ' num2str(error_SDP)]);
    tol = abs(SDP_obj_all(counter+1) - SDP_obj_all(counter)); % update obj tol
    str = ['tol = ' num2str(tol)];
   disp(str);
    
    soft_output(1:unlabelednum,counter+1)=x_valid(labelednum+1:end,:);
    pause(0.00001);
   
end
ck=c_k;

class_SDP_temp_binary = sign(x_valid);
class_SDP_temp_binary(initial_label_index) = [];

class_SDP_temp = zeros(size(class_SDP_temp_binary,1),size(class_SDP_temp_binary,2));

%% ======FINAL LABELS AFTER FEATURE GRAPH UPDATE======
class_SDP_temp(class_SDP_temp_binary==1) = 1;
class_SDP_temp(class_SDP_temp_binary==-1) = -1;
%% ===================================================

%%=========================================================================
end