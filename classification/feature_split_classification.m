clear
clc

[labeled_event,~]=xlsread('feature table');
classifier_model=['SVM'];

total_run=50;
score_GTV=cell(total_run,1);
score_normGLR=cell(total_run,1);
confusion_matrix_GTV=cell(total_run,1);
confusion_matrix_normGLR=cell(total_run,1);
sensitive_GTV=zeros(total_run,4);
sensitive_normGLR=zeros(total_run,4);

for i_run=1:total_run

P = 0.70 ;

[r_class1_index,c_class1_index]=find(labeled_event(:,7)==1);
[r_class2_index,c_class2_index]=find(labeled_event(:,7)==2);
[r_class3_index,c_class3_index]=find(labeled_event(:,7)==3);
[r_class4_index,c_class4_index]=find(labeled_event(:,7)==4);

class1=labeled_event(r_class1_index,:);
class2=labeled_event(r_class2_index,:);
class3=labeled_event(r_class3_index,:);
class4=labeled_event(r_class4_index,:);

[m,n] = size(class1) ;
idx = randperm(m)  ;
Training_class1 = class1(idx(1:round(P*m)),:) ; 
Testing_class1 = class1(idx(round(P*m)+1:end),:) ;

[m,n] = size(class2) ;
idx = randperm(m)  ;
Training_class2 = class2(idx(1:round(P*m)),:) ; 
Testing_class2 = class2(idx(round(P*m)+1:end),:) ;

[m,n] = size(class3) ;
idx = randperm(m)  ;
Training_class3 = class3(idx(1:round(P*m)),:) ; 
Testing_class3 = class3(idx(round(P*m)+1:end),:) ;

[m,n] = size(class4) ;
idx = randperm(m)  ;
Training_class4 = class4(idx(1:round(P*m)),:) ; 
Testing_class4 = class4(idx(round(P*m)+1:end),:) ;

Training=[Training_class1;Training_class2;Training_class3;Training_class4];
Testing=[Testing_class1;Testing_class2;Testing_class3;Testing_class4];

train_feature=Training(:,8:end);
train_label=Training(:,7);
test_feature=Testing(:,8:end);
test_label=Testing(:,7);


%% inverse q folder cross validation RF SVM GLR norm-GLR GTV 
final_batch_classification_result=zeros(length(test_label),4);

    
    rescale_boxcox_train_feature = zeros(size(train_feature));
    rescale_boxcox_test_feature = zeros(size(test_feature));
    rescale_boxcox_train_feature_norm = zeros(size(train_feature));
    rescale_boxcox_test_feature_norm=zeros(size(test_feature));
    min_rescale_boxcox_train_feature = zeros(121,1);
    max_rescale_boxcox_train_feature = zeros(121,1);
            
for iii=1:119
      rescale_boxcox_train_feature(:,iii)=train_feature(:,iii);
        rescale_boxcox_test_feature(:,iii)=test_feature(:,iii);

min_rescale_boxcox_train_feature(iii)=min(rescale_boxcox_train_feature(:,iii));
max_rescale_boxcox_train_feature(iii)=max(rescale_boxcox_train_feature(:,iii)-min_rescale_boxcox_train_feature(iii));
rescale_boxcox_train_feature_norm(:,iii)=(rescale_boxcox_train_feature(:,iii)-min_rescale_boxcox_train_feature(iii))/max_rescale_boxcox_train_feature(iii); % (max_rescale_boxcox_train_feature(iii)-min_rescale_boxcox_train_feature(iii))
rescale_boxcox_test_feature_norm(:,iii)=(rescale_boxcox_test_feature(:,iii)-min_rescale_boxcox_train_feature(iii))/max_rescale_boxcox_train_feature(iii); % (max_rescale_boxcox_train_feature(iii)-min_rescale_boxcox_train_feature(iii))
end
    
    %% PCA feature extraction begin
     train_feature_classifier_PCA=[];
     test_feature_classifier_PCA=[];

        train_feature_classifier_PCA=rescale_boxcox_train_feature_norm;
        test_feature_classifier_PCA=rescale_boxcox_test_feature_norm;

%     %% PCA 
    train_feature_PCA=[];
    test_feature_PCA=[];
    features=train_feature_classifier_PCA;
    [i1,j]=size(features);
    [coeff,score,latent,tsquared,explained,mu1] = pca(features);
    pca_own.coeff=coeff;
    pca_own.score=score;
    pca_own.latent=latent;
    pca_own.tsquared=tsquared;
    pca_own.explained=explained;
    pca_own.mu1=mu1;
    
    train_feature_PCA = pca_own.score;
    [r_test_feature,c_test_feature]=size(test_feature_classifier_PCA);
    test_feature_PCA=(test_feature_classifier_PCA-repmat(pca_own.mu1,r_test_feature,1))*pca_own.coeff;
    
    B_pca = cumsum(pca_own.explained);
    cum_pca = B_pca/sum(pca_own.explained);
    diff_cum_pca = abs(cum_pca-0.95);
    [pca_feature_num, pca_value]=find(diff_cum_pca==min(diff_cum_pca));
    
    %% boxcox norm std PCA
    
    rescale_boxcox_train_feature_PCA = zeros(size(train_feature_PCA));
    rescale_boxcox_test_feature_PCA = zeros(size(test_feature_PCA));
    rescale_boxcox_train_feature_norm_PCA = zeros(size(train_feature_PCA));
    rescale_boxcox_test_feature_norm_PCA=zeros(size(test_feature_PCA));
    rescale_boxcox_train_feature_std_PCA = zeros(size(train_feature_PCA));
    rescale_boxcox_test_feature_std_PCA=zeros(size(test_feature_PCA));
    min_rescale_boxcox_train_feature_PCA = zeros(121,1);
    max_rescale_boxcox_train_feature_PCA = zeros(121,1);
    std_rescale_boxcox_train_feature_PCA = zeros(121,1);
    
rescale_boxcox_train_feature_PCA=train_feature_PCA;
rescale_boxcox_test_feature_PCA=test_feature_PCA;

    for iiij=1:119
              
min_rescale_boxcox_train_feature_PCA(iiij)=min(rescale_boxcox_train_feature_PCA(:,iiij));
max_rescale_boxcox_train_feature_PCA(iiij)=max(rescale_boxcox_train_feature_PCA(:,iiij)-min_rescale_boxcox_train_feature_PCA(iiij));
rescale_boxcox_train_feature_norm_PCA(:,iiij)=(rescale_boxcox_train_feature_PCA(:,iiij)-min_rescale_boxcox_train_feature_PCA(iiij))/max_rescale_boxcox_train_feature_PCA(iiij);
rescale_boxcox_test_feature_norm_PCA(:,iiij)=(rescale_boxcox_test_feature_PCA(:,iiij)-min_rescale_boxcox_train_feature_PCA(iiij))/max_rescale_boxcox_train_feature_PCA(iiij);
    end
    
        train_feature_classifier=rescale_boxcox_train_feature_norm_PCA(:,1:pca_feature_num(1));
        test_feature_classifier=rescale_boxcox_test_feature_norm_PCA(:,1:pca_feature_num(1));

    
    num_top_coeff=pca_feature_num(1);
    C=[];sen_quake=[];sen_seismic=[];sen_rockfall=[];sen_noise=[];
    [final_score_output]=final_output_feature(train_feature_classifier,test_feature_classifier,train_label,test_label,classifier_model,num_top_coeff);
    
    [~,pred_label]=max(final_score_output,[],2);
    C=confusionmat(test_label, pred_label)
    sen_quake=C(1,1)/sum(C(1,:))
    sen_seismic=C(2,2)/sum(C(2,:))
    sen_rockfall=C(3,3)/sum(C(3,:))
    sen_noise=C(4,4)/sum(C(4,:))
    sensitive_C_PCA_SVM(i_run,1)=sen_quake;
    sensitive_C_PCA_SVM(i_run,2)=sen_seismic;
    sensitive_C_PCA_SVM(i_run,3)=sen_rockfall;
    sensitive_C_PCA_SVM(i_run,4)=sen_noise;
    confusion_matrix_PCA_SVM{i_run}=C;
    score_PCA_SVM{i_run}=final_score_output;
    

    
end
save_confusion_PCA_SVM=['classifier_PCA_SVM_withoutgraphbf_30test_confusion_matrix.mat'];
save_sensitive_PCA_SVM=['classifier_PCA_SVM_cvx_withoutgraphbf_30test_sensitive.mat'];
save_score_PCA_SVM=['score_PCA_SVM_withoutgraphbf_30test.mat'];

save(save_confusion_PCA_SVM,'confusion_matrix_PCA_SVM');
save(save_sensitive_PCA_SVM,'sensitive_C_PCA_SVM');
save(save_score_PCA_SVM,'score_PCA_SVM');


