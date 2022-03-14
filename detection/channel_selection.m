
clear
clc
close all

%% window segment

[detected_stack_signal,~]=xlsread('detection result.xlsx','detected MCM dataset1');


[a,b]=find(whole_detected_dataset3(:,5)==2);
whole_detected_dataset3(a,5)=whole_detected_dataset3(a,5)+1;

dataset_event=[num2str(detected_stack_signal(:,1)),num2str(detected_stack_signal(:,2),'%02d'),num2str(detected_stack_signal(:,3),'%02d'),num2str(detected_stack_signal(:,4),'%02d')];
dataset_event = double(string(dataset_event))*1000000+(detected_stack_signal(:,6)+detected_stack_signal(:,5)*60)*250;
[r_w_data,~]=size(detected_stack_signal);

dataset_name=dataset_event;
path=['\...'];
save_path=['\...'];

event_start_index_vector=zeros(r_w_data,1);

for i=1:r_w_data
   
    segment=[];
    file_name=['bfstart_hour_',num2str(dataset_name(event_start_index_vector(i),1)),'-',num2str(dataset_name(event_start_index_vector(i),2),'%02d'),'-',num2str(dataset_name(event_start_index_vector(i),3),'%02d'),'-',num2str(dataset_name(event_start_index_vector(i),4),'%02d'),'.mat'];
    start_index=(dataset_name(event_start_index_vector(i),5)*60+dataset_name(event_start_index_vector(i),6))*250;
    
    load([path,file_name])
    save_file_name=['window_',num2str(dataset_name(event_start_index_vector(i),1)),'-',num2str(dataset_name(event_start_index_vector(i),2),'%02d'),'-',num2str(dataset_name(event_start_index_vector(i),3),'%02d'),'-',num2str(dataset_name(event_start_index_vector(i),4),'%02d'),'-',num2str(dataset_name(event_start_index_vector(i),5),'%02d'),'-',num2str(dataset_name(event_start_index_vector(i),6),'%02d'),'.mat'];
    
    [r_size,c_size]=size(norm_data_BF);
    
channel_index=1:c_size;

    if start_index+2.5*60*250>900001
        segment=norm_data_BF(start_index-2.5*60*250:900001,channel_index);
        segment = [segment;zeros(fix(start_index+2.5*60*250-900001),length(channel_index))];
    elseif  start_index-2.5*60*250<=0  
        segment=norm_data_BF(1:start_index+2.5*60*250,channel_index);
        segment = [zeros(fix(abs(start_index-2.5*60*250)),length(channel_index));segment];
    else
        segment=norm_data_BF(start_index-2.5*60*250:start_index+2.5*60*250,channel_index);
    end
    save([save_path,save_file_name],'segment')
    
end




%%
dataset_path=[save_path];
d = dir(dataset_path);
window_length=100;
selected_channel=0;
final_file=0;
for i=3:length(d)
    disp(i)
    ratio = zeros(1,6);
 file_name = d(i).name;
   load([dataset_path,file_name])
    year = double(string(file_name(8:11)));
    month = double(string(file_name(13:14)));
    day = double(string(file_name(16:17)));
    hour = double(string(file_name(19:20)));
    min = double(string(file_name(22:23)));
    second = double(string(file_name(25:27)));
    
    [time_length,channel_length]=size(segment);
    
    % replace zero
    for jj=1:channel_length
        se=zeros(time_length,1);
        sub_segment=segment(:,jj);
        sub_segment(find(abs(sub_segment<10^(-15))))=mean(sub_segment);
        se=sub_segment;
        segment(1:time_length,jj)=(se-mean(se))/std(se);
    end
    
    
    power_window=0;
    for j=1:fix(time_length/window_length)
        for channel_index = 1:channel_length
        power_window(j,channel_index) = sum(segment(1+window_length*(j-1):j*window_length,channel_index).^2);
        end
    end
    
    for channel=1:channel_length
    [f,x] = ecdf(power_window(:,channel));
[row_index,col_index] = find(f>0.0001&f<0.9999);
ratio(channel) = mean(x(row_index(end-5:end)))/mean(x(row_index(1:6)));
    s=1;
    end

    [r_index,c_index]=find(ratio==max(ratio));
    selected_channel(i-2)=c_index;
    final_file(i-2,1:7) = [year,month,day,hour,min,second,c_index];
end

% c_index is the selected channel index


