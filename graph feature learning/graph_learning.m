%% 
clear;
clc;
close all;


class=[1 2 3 4];
final_score=[];
final_ck=zeros(4,119);
global labelednum unlabelednum 

% signal split
[labeled_event,~]=xlsread('feature table');
classifier_model=['GLR_CVX'];

total_run=50;
confusion_matrix=cell(total_run,1);
sensitive=zeros(total_run,4);

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


for class_index=1:4
    train_feature=[];
    train_label=[];
    test_feature=[];
    test_label=[];
    
    train_feature=Training(:,8:end);
train_label=Training(:,7);
test_feature=Testing(:,8:end);
test_label=Testing(:,7);
    
if class_index==1
train_label(train_label==1)=1;
train_label(train_label==2)=-1;
train_label(train_label==3)=-1;
train_label(train_label==4)=-1;
test_label(test_label==1)=1;
test_label(test_label==2)=-1;
test_label(test_label==3)=-1;
test_label(test_label==4)=-1;
end

if class_index==2
train_label(train_label==1)=-1;
train_label(train_label==2)=1;
train_label(train_label==3)=-1;
train_label(train_label==4)=-1;
test_label(test_label==1)=-1;
test_label(test_label==2)=1;
test_label(test_label==3)=-1;
test_label(test_label==4)=-1;
end

if class_index==3
train_label(train_label==1)=-1;
train_label(train_label==2)=-1;
train_label(train_label==3)=1;
train_label(train_label==4)=-1;
test_label(test_label==1)=-1;
test_label(test_label==2)=-1;
test_label(test_label==3)=1;
test_label(test_label==4)=-1;
end

if class_index==4
train_label(train_label==1)=-1;
train_label(train_label==2)=-1;
train_label(train_label==3)=-1;
train_label(train_label==4)=1;
test_label(test_label==1)=-1;
test_label(test_label==2)=-1;
test_label(test_label==3)=-1;
test_label(test_label==4)=1;
end

  

    rescale_boxcox_train_feature_norm = zeros(size(train_feature));
    rescale_boxcox_test_feature_norm=zeros(size(test_feature));
    min_rescale_boxcox_train_feature = zeros(119,1);
    max_rescale_boxcox_train_feature = zeros(119,1);
    

    for iii=1:119

min_rescale_boxcox_train_feature(iii)=min(train_feature(:,iii));
max_rescale_boxcox_train_feature(iii)=max(train_feature(:,iii)-min_rescale_boxcox_train_feature(iii));
rescale_boxcox_train_feature_norm(:,iii)=(train_feature(:,iii)-min_rescale_boxcox_train_feature(iii))/max_rescale_boxcox_train_feature(iii);
rescale_boxcox_test_feature_norm(:,iii)=(test_feature(:,iii)-min_rescale_boxcox_train_feature(iii))/max_rescale_boxcox_train_feature(iii);
    end
  

  train_feature_classifier=rescale_boxcox_train_feature_norm;
    test_feature_classifier=rescale_boxcox_test_feature_norm;
    
    class_test = test_label;
    labelednum = length(train_label);
    unlabelednum = length(test_label);
    feature_train_test =[train_feature_classifier;test_feature_classifier];
    initial_label=[train_label;zeros(unlabelednum,1)];
    initial_label_index = logical([ones(labelednum,1); zeros(unlabelednum, 1)]);
    class_train_test=[train_label;test_label];
    
    
sprintf(['the classification result of class',num2str(class_index)])

%% graph learning setting

disp(unlabelednum)
disp('=========start========');
num_feature=119;
S_lower=0;
S_upper=5000;
c_k_initially_set =1190;
tol_set=10^(-4);
tol_set_pg=10^(-4);

initial_label=[train_label;zeros(unlabelednum,1 )];
initial_label_index = logical([ones(labelednum,1); zeros(unlabelednum, 1)]);
count_class0_partA=length(find(test_label==-1));

alpha=1;
x_known=[train_label;ones(length(class_train_test)-length(train_label),1)*0];

sigma=1;

[ x_valid, class_SDP_temp, SDP_error, ck_A, whole_final_ck] = sdp_binary_GU_oao_L_constant_sign_o_norm_replace_GLR_GTV_L_norm( class_test, ...
    c_k_initially_set, tol_set, tol_set_pg, S_lower, S_upper, ...
    feature_train_test, initial_label, initial_label_index, class_train_test,alpha,train_label,sigma);
save(['normGLR_ck_init1_class',num2str(class_index),'_ck.mat'],'whole_final_ck')
save(['normGLR_ck_init1_class',num2str(class_index),'_label.mat'],'x_valid')


end


diary off