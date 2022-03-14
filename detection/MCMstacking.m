clear
clc
close all


path=['\...']; % dataset stored path
d=dir(path);
save_path=['\...'];


for d_index=3:length(d)
      disp(d_index)
data=load([path,d(d_index).name]);

start_time=1;
[end_time,channel_num]=size(data.norm_data_BF);


norm_data=zeros(end_time,channel_num);
norm_data_BF=zeros(end_time,channel_num);
window_length=50; % window length selection is a trade off of the waveform resolution and time spend
fc1=5;
fc2 =100;
fs = 250;
[b,a] = butter(8,[fc1 fc2]/(fs/2));
for i=1:channel_num
    norm_data(:,i)=data.norm_data_BF(:,i)-mean(data.norm_data_BF(:,i)); % set mean value as zero
    norm_data_BF(:,i) = filter(b,a,norm_data(:,i)); % filteing
end


time_windows=zeros(end_time,channel_num);
MCM=zeros(end_time-window_length,channel_num);

for i =1:window_length:(end_time-window_length)
    
    for j =1:channel_num
    time_windows(i:i+window_length,j)=norm_data_BF(i:i+window_length,j)-mean(norm_data_BF(i:i+window_length,j));
    end
    index=0;
for j=1:channel_num-2
    for num=j+1:channel_num-1
        for m=num+1:channel_num
        index=index+1;
a=j;
b=num;
c=m;
upper=time_windows(i:i+window_length,a).*time_windows(i:i+window_length,b).*time_windows(i:i+window_length,c);
MCM(i,index)=sum(abs(upper));
        end
    end
end
end
MCM_stacking=(sum(MCM,2));

save([save_path,'MCM',d(d_index).name],'MCM_stacking')
end


