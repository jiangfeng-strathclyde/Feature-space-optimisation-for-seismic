function [final_score_output]=final_output_feature(train_data,test_data,train_label,test_label,classifier_model,num_top_coeff)

X =train_data;
z =test_data;
class=[1 2 3 4];
[row_test,col_test]=size(z);


final_classification=zeros(row_test,1);
for num_class=1:4
X_train=[];
Y_train=[];
X_test=[];
Y_test=[];
[row_train,~]=size(X);
Y_train=train_label;
[aa,bb]=find(Y_train~=class(num_class));
Y_train(aa)=-1;
[aaa,bba]=find(Y_train==class(num_class));
Y_train(aaa)=1;
X_train = X;
[row_test,~]=size(z);
Y_test=test_label;
[aa,bb]=find(Y_test~=class(num_class));
Y_test(aa)=-1;
[aaa,bba]=find(Y_test==class(num_class));
Y_test(aaa)=1;
X_test=z;

switch classifier_model
    case 'SVM'        
          SVMModel = fitcsvm(X_train(:,1:num_top_coeff),Y_train,'KernelFunction','rbf','KernelScale','auto');
        [predict_label,Scores] = predict(SVMModel,X_test(:,1:num_top_coeff));
        [row_score,~]=size(Scores);
        for i=1:row_score
            if Scores(i,1)>=Scores(i,2)
                final_score(i)=-Scores(i,1);
            else if Scores(i,1)<Scores(i,2)
        final_score(i)=Scores(i,2);
                end
            end
        end
       final_score=final_score;
        
    case 'RF'
        Factor = TreeBagger(50, X_train(:,1:num_top_coeff),Y_train);
        [Predicted_label,Scores] = predict(Factor, X_test(:,1:num_top_coeff));
        predicted=zeros(size(Predicted_label,1),size(Predicted_label,2));
        predicted=str2double(Predicted_label);
        [row_score,~]=size(Scores);
        for i=1:row_score
            
            if Scores(i,1)>=Scores(i,2)
                final_score(i)=-Scores(i,1);
            else if Scores(i,1)<Scores(i,2)
        final_score(i)=Scores(i,2);
                end
            end
        end
        final_score=final_score;
        
case 'GLR_CVX'
    
    feature_number=num_top_coeff;
        S_label=[Y_train;Y_test];    
        num_train_class1 = length(find(Y_train==1));
        num_train_class0=length(Y_train)-num_train_class1;
        training_features=X_train(:,1:feature_number);
        testing_features=X_test(:,1:feature_number);
        [i1,j1]=size(training_features);
        [i2,j2]=size(testing_features);
        total_features=zeros(i1+i2,feature_number);
        total_features=[training_features;testing_features];
        s=zeros(i1+i2);
        A=zeros(i1+i2);
        As=zeros(i1+i2);
        ss=cell(feature_number,1); 
        FD=cell(feature_number,1);
        diff=[];
        
       for num_features=1:feature_number

             diff{num_features}=(total_features(:,num_features)-total_features(:,num_features).').^2;   % adjacency matrix of each segment
            
             %% sigma setup

            feature_diff=diff{num_features}*(1);

            W_class1_same= feature_diff(1:num_train_class1,1:num_train_class1);
            W_class1_diff= feature_diff(1:num_train_class1,num_train_class1+1:(num_train_class1+num_train_class0));
            W_class0_same= feature_diff(num_train_class1+1:(num_train_class1+num_train_class0),num_train_class1+1:(num_train_class1+num_train_class0));

            W_class1_same=triu(W_class1_same); 
            W_class1_same=W_class1_same+diag(-diag(W_class1_same));
            W_class1_same=W_class1_same(:);
            W_class1_same(W_class1_same==0)=[];

            W_class0_same=triu(W_class0_same); 
            W_class0_same=W_class0_same+diag(-diag(W_class0_same));
            W_class0_same=W_class0_same(:);
            W_class0_same(W_class0_same==0)=[];
            W_class1_diff=W_class1_diff(:);


            same_label=[W_class0_same;W_class1_same];
            diff_label=W_class1_diff;


            A=same_label';
            B=diff_label';
            n=length(A);
            m=length(B);

            syms p
            p=@(x) (A*exp((-A')./x.^2)/n - B*exp((-B')./x.^2)/m);
            [y]=erfenfa_new_use2(p,0.1,1.2,0.001,a);
            sigma(num_features)=y; 

            s=-diff{num_features}/(y)^2;   % adjacency matrix of each segment
            ss{num_features}=s;
            As=As+s;
       end
       
        half_D=[];
        W=exp(As);
        W(W == diag(W)) = 0;    
        D=diag(sum(W));%     diagonal degree matrix
        % remove zeros values from D outliers
        D_diag=diag(D);
        [z_r,z_c]=find(D_diag==0);
        D(z_r,z_r)=mean(D_diag);
           D_diag=diag(D);     
        L=D-W;
        L(isnan(L))=0;
        L(isinf(L))=0;
                    
% CVX learned label
                            n=length(S_label);
                            cvx_begin sdp
                                variable nu(n)
                                maximize ( -sum(nu) )
                                L+ diag(nu) >= 0;
                            cvx_end
                            
                            new_L=L+diag(nu);
                            W3=new_L(i1+1:end,i1+1:end);
                            W1=new_L(1:i1,1:i1);
                            W2=new_L(1:i1,i1+1:end);
                            x1=S_label(1:i1);
                            x2=S_label(i1+1:end);

                            B=W2'*x1;
                            C=x1'*W2;
                            new_x=-inv(W3)*(B+C')/2;
                            
                            final_score=new_x;
                        
                      
    case 'normGLR_CVX'
        
       
        feature_number=num_top_coeff;
        S_label=[Y_train;Y_test];    
        num_train_class1 = length(find(Y_train==1));
        num_train_class0=length(Y_train)-num_train_class1;
        training_features=X_train(:,1:feature_number);
        testing_features=X_test(:,1:feature_number);
        [i1,j1]=size(training_features);
        [i2,j2]=size(testing_features);
        total_features=zeros(i1+i2,feature_number);
        total_features=[training_features;testing_features];
        s=zeros(i1+i2);
        A=zeros(i1+i2);
        As=zeros(i1+i2);
        ss=cell(feature_number,1); 
        FD=cell(feature_number,1);
        diff=[];
        
       for num_features=1:feature_number

             diff{num_features}=(total_features(:,num_features)-total_features(:,num_features).').^2;   % adjacency matrix of each segment

             %% sigma setup

            feature_diff=diff{num_features}*(1);

            W_class1_same= feature_diff(1:num_train_class1,1:num_train_class1);
            W_class1_diff= feature_diff(1:num_train_class1,num_train_class1+1:(num_train_class1+num_train_class0));
            W_class0_same= feature_diff(num_train_class1+1:(num_train_class1+num_train_class0),num_train_class1+1:(num_train_class1+num_train_class0));

            W_class1_same=triu(W_class1_same); 
            W_class1_same=W_class1_same+diag(-diag(W_class1_same));
            W_class1_same=W_class1_same(:);
            W_class1_same(W_class1_same==0)=[];

            W_class0_same=triu(W_class0_same); 
            W_class0_same=W_class0_same+diag(-diag(W_class0_same));
            W_class0_same=W_class0_same(:);
            W_class0_same(W_class0_same==0)=[];
            W_class1_diff=W_class1_diff(:);


            same_label=[W_class0_same;W_class1_same];
            diff_label=W_class1_diff;


            A=same_label';
            B=diff_label';
            n=length(A);
            m=length(B);

            syms p
            p=@(x) (A*exp((-A')./x.^2)/n - B*exp((-B')./x.^2)/m);
            [y]=erfenfa_new_use2(p,0.1,1.2,0.001,a);
            sigma(num_features)=y; 

            s=-diff{num_features}/(y)^2;   % adjacency matrix of each segment
            ss{num_features}=s;
            As=As+s;
       end
       
               half_D=[];
        W=exp(As);
        W(W == diag(W)) = 0;    
        D=diag(sum(W));%     diagonal degree matrix
        % remove zeros values from D outliers
        D_diag=diag(D);
        [z_r,z_c]=find(D_diag==0);
        D(z_r,z_r)=mean(D_diag);
           D_diag=diag(D);     
        L=D-W;
        L(isnan(L))=0;
        L(isinf(L))=0;

        half_D=zeros(size(D));
        half_D(logical(eye(size(half_D))))=D_diag.^(-1/2);
        L_norm=half_D*L*half_D;
            T=isnan(L_norm);
            [row_L_norm,~]=size(L_norm);
                if sum(sum(T))~=0
                    L_norm=ones(row_L_norm);
                end
        TF_normGLR=isinf(L_norm(i1+1:(i1+i2),i1+1:(i1+i2)));       
        [row_inf,col_inf]=find(L_norm==Inf);
        if isempty(row_inf)~=1
        L_norm(row_inf,col_inf)=1;
        end
        [row_nan,col_nan]=find(isnan(L_norm));
        if isempty(row_nan)~=1
        L_norm(row_nan,col_nan)=0;
        end
          L_norm(isnan(L_norm))=0;
        L_norm(isinf(L_norm))=0;
        
n=[];new_L=[];W1=[];W2=[];W3=[];x1=[];x2=[];B=[];C=[];new_x=[];
                            n=length(S_label);
                            cvx_begin sdp
                                variable nu(n)
                                maximize ( -sum(nu) )
                                L_norm+ diag(nu) >= 0;
                            cvx_end
                            
                            new_L=L_norm+diag(nu);
                            W3=new_L(i1+1:end,i1+1:end);
                            W1=new_L(1:i1,1:i1);
                            W2=new_L(1:i1,i1+1:end);
                            x1=S_label(1:i1);
                            x2=S_label(i1+1:end);

                            B=W2'*x1;
                            C=x1'*W2;
                            new_x=-inv(W3)*(B+C')/2;
                            
                            final_score=new_x;
                       

%% GTV_CVX

    case 'GTV_CVX' 
 feature_number=num_top_coeff;
        S_label=[Y_train;Y_test];    
       num_train_class1 = length(find(Y_train==1));
        num_train_class0=length(Y_train)-num_train_class1;
        training_features=X_train(:,1:feature_number);
        testing_features=X_test(:,1:feature_number);
        [i1,j1]=size(training_features);
        [i2,j2]=size(testing_features);
        total_features=zeros(i1+i2,feature_number);
        total_features=[training_features;testing_features];
        s=zeros(i1+i2);
        A=zeros(i1+i2);
        As=zeros(i1+i2);
        ss=cell(feature_number,1); 
        FD=cell(feature_number,1);
        diff=[];
        
       for num_features=1:feature_number

             diff{num_features}=(total_features(:,num_features)-total_features(:,num_features).').^2;   % adjacency matrix of each segment
            
             %% sigma setup

            feature_diff=diff{num_features}*(1);

            W_class1_same= feature_diff(1:num_train_class1,1:num_train_class1);
            W_class1_diff= feature_diff(1:num_train_class1,num_train_class1+1:(num_train_class1+num_train_class0));
            W_class0_same= feature_diff(num_train_class1+1:(num_train_class1+num_train_class0),num_train_class1+1:(num_train_class1+num_train_class0));

            W_class1_same=triu(W_class1_same); 
            W_class1_same=W_class1_same+diag(-diag(W_class1_same));
            W_class1_same=W_class1_same(:);
            W_class1_same(W_class1_same==0)=[];

            W_class0_same=triu(W_class0_same); 
            W_class0_same=W_class0_same+diag(-diag(W_class0_same));
            W_class0_same=W_class0_same(:);
            W_class0_same(W_class0_same==0)=[];
            W_class1_diff=W_class1_diff(:);


            same_label=[W_class0_same;W_class1_same];
            diff_label=W_class1_diff;


            A=same_label';
            B=diff_label';
            n=length(A);
            m=length(B);

            syms p
            p=@(x) (A*exp((-A')./x.^2)/n - B*exp((-B')./x.^2)/m);
            [y]=erfenfa_new_use2(p,0.1,1.2,0.001,a);
            sigma(num_features)=y; 

            s=-diff{num_features}/(y)^2;   % adjacency matrix of each segment
            ss{num_features}=s;
            As=As+s;
       end
       
        half_D=[];
        W=exp(As);
        W(W == diag(W)) = 0;    
        D=diag(sum(W));%     diagonal degree matrix
        % remove zeros values from D outliers
        D_diag=diag(D);
        [z_r,z_c]=find(D_diag==0);
        D(z_r,z_r)=mean(D_diag);
           D_diag=diag(D);     
        L=D-W;
        L(isnan(L))=0;
        L(isinf(L))=0;
  



I=[];A=[]; n=[];new_L=[];W1=[];W2=[];W3=[];x1=[];x2=[];B=[];C=[];new_x=[];
                             I =eye(i1+i2);
                            A=(I-W)'*(I-W);
                           A(isnan(A))=0;
                           A(isinf(A))=0;
                        
                            n=length(S_label);
                            cvx_begin sdp
                                variable nu(n)
                                maximize ( -sum(nu) )
                                A+ diag(nu) >= 0;
                            cvx_end
                            
                            new_L=A+diag(nu);
                            W3=new_L(i1+1:end,i1+1:end);
                            W1=new_L(1:i1,1:i1);
                            W2=new_L(1:i1,i1+1:end);
                            x1=S_label(1:i1);
                            x2=S_label(i1+1:end);

                            B=W2'*x1;
                            C=x1'*W2;
                            new_x=-inv(W3)*(B+C')/2;
                            
                            final_score=new_x;

   
    otherwise

        disp('other value')
end
 final_classification(:,num_class)=final_score;
end
final_score_output=final_classification;
end
        
      
