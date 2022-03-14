%% for detected true event stacked_channel
clear
clc
close all

dataset1_path=['...\dataset1\'];
dataset2_path=['...\dataset2\'];
dataset3_path=['...\dataset3\'];

a_segment_length=500;
a_segment_num=50;
level_count=1;
event_interval=50;
window_length = 50;


for dataset_index=1:3
    level_count=1;
    if dataset_index==1
        save_file=['siteA_dataset1.mat'];
        dataset_path=dataset1_path;
    end
        if dataset_index==2
        save_file=['siteA_dataset2.mat'];
        dataset_path=dataset2_path;
        end
        if dataset_index==3
        save_file=['siteA_dataset3.mat'];
        dataset_path=dataset3_path;
        end
    
    final_result=zeros(1,9);
d=dir(dataset_path);
for signal_index=3:length(d)
     disp(signal_index)
   pairs_channel1=zeros(1,2);
   diff=0;
    load([dataset_path,d(signal_index).name])
    date_name=d(signal_index).name(17:end-4);
    year=double(string(date_name(1:4)));
    month=double(string(date_name(6:7)));
    day=double(string(date_name(9:10)));
    hour=double(string(date_name(12:13)));
    
    a=MCM_stacking;
    a(a==0)=mean(abs(a));
    a=abs(a);
    r_duration=[];
    count_index=0;

%% threshold setup
while isempty(r_duration)==1

randon_a_segment_start_index=randi([1 length(a)-a_segment_length], a_segment_num,1);  
      
for ii=1:a_segment_num
    a_segment=a(randon_a_segment_start_index(ii):randon_a_segment_start_index(ii)+a_segment_length);
    [D12 PD12] = allfitdist(a_segment,'PDF');
    nu_value(dataset_index,ii)=PD12{1,1}.nu;
    
    nu= nu_value(dataset_index,ii);
    xbar=mean(a_segment);
    se=std(a_segment);
    crit = tinv(0.99999,nu);
    ci = xbar + crit*se/sqrt(a_segment_length);

    lamda(dataset_index,ii)=ci;
      end
  
      [row,col]=size(lamda);
new_shape = reshape(lamda,[1,row*col]);
[f,x] = ecdf(new_shape);
[row_index,col_index] = find(f>0.1&f<0.9);
lamda_value = max(x(row_index));
the = lamda_value;

%% detection 
    [row1,~]=find(a>the);
    picks=row1;  
    
    if length(picks)>1
    A=picks;
    c1 = 1;
arrset = cell(0,0);
while(c1<numel(A))
    c2 = 0;
    while (c1+c2+1<=numel(A)&&A(c1)+c2+1==A(c1+c2+1))
        c2 = c2+1;
    end
    if(c2>=0)
        arrset= [arrset;(A(c1:1:c1+c2))];
    end
    c1 = c1 + c2 +1;
end


if isempty(arrset)~=1
    if length(arrset)==1
        pairs_channel1(1,1)=arrset{1,1}(1);
        pairs_channel1(1,2)=arrset{1,1}(end);
    else 
        event_index=0;
        diff=zeros(1,length(arrset)-1);
        for jj=2:length(arrset)
            diff(jj-1)=arrset{jj,1}(1)-arrset{jj-1,1}(end);
        end
        
        [r_event_interval,c_event_interval]=find(diff>event_interval);
        event_index=[0,c_event_interval,length(arrset)];
        
        if isempty(c_event_interval)==1
        pairs_channel1(1,1)=arrset{1,1}(1);
        pairs_channel1(1,2)=arrset{end,1}(end);
        else 
            for i_diff = 1:length(c_event_interval)+1
                pairs_channel1(i_diff,1)=arrset{event_index(i_diff)+1,1}(1);
                pairs_channel1(i_diff,2)=arrset{event_index(i_diff+1),1}(end);
            end
        end
    
end
   
end

old_duration=pairs_channel1(:,2)-pairs_channel1(:,1);
pairs_channel1(find(old_duration==0),:)=[];
duration=pairs_channel1(:,2)-pairs_channel1(:,1);
[r_duration,c_duration]=find(duration>100);
    end
count_index=count_index+1;

if count_index>30
    pairs_channel1=[0,0];
    break
end

end
[row_pairs,~]=size(pairs_channel1);
final_result(level_count:level_count+row_pairs-1,1)=year;
final_result(level_count:level_count+row_pairs-1,2)=month;
final_result(level_count:level_count+row_pairs-1,3)=day;
final_result(level_count:level_count+row_pairs-1,4)=hour;
final_result(level_count:level_count+row_pairs-1,5)=floor(pairs_channel1(:,1)/250/60);
final_result(level_count:level_count+row_pairs-1,6)=(pairs_channel1(:,1)/250 - floor(pairs_channel1(:,1)/250/60)*60);
final_result(level_count:level_count+row_pairs-1,9)=pairs_channel1(:,2)-pairs_channel1(:,1);
final_result(level_count:level_count+row_pairs-1,7)=pairs_channel1(:,1);
final_result(level_count:level_count+row_pairs-1,8)=pairs_channel1(:,2);
level_count=level_count+row_pairs;

end

save_path=['...\detectionwith6channel\'];
save([save_path,save_file],'final_result')

end
