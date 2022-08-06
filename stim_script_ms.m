
%% post-stimuli analysis NETWORK RESPONSE
files = dir('*_cSpikes_L0_RP2*.mat')
for i=1:length(files); %for the list of files
    fileID    =   files(i).name;  % at any one time the file in memory 
   %is selected from the list (i)
   data =  load(fileID);
try  
spikes= data.cSpikes; %or cSpikes
catch
    spikes = data.spikes
end
channels = data.channels;
recordDuration = round(length(spikes)); %in samples
 %downSpikeMatrix = downSampleSum_org(spikeMatrix, recordDuration * 1/fs);
  % fs = 25; %mHz
   %duration_ms = 300000; %ms
%% truncate data into 1 ms bins
   duration_ms = 300000; 
   %downFactor = 25000; %per s
   downFactor_ms = 25;
   electrode_no = 60;
   %downSpikeMatrix = downSampleSum(full(spikes), recordDuration/downFactor);
   %in s
   downTrain = reshape(full(spikes), [], duration_ms, electrode_no);
    downTrain = sum(downTrain, 1); 
    outTrain = reshape(downTrain, duration_ms, electrode_no);

index = 1:60;
ntrials = 60; %number of trials  per te 5 minute recording = 60
electrode_no = 60;

%% get matrix of only post-stim times (100ms)
Post_Stim_Matrix = zeros(ntrials*100, electrode_no);
skip_durationms = 5000 %5000/ 10 ms to skip to next 5 s ( when is pulse administered)
newpairB = 0:100:length(Post_Stim_Matrix);
skipA = 0:skip_durationms:length(outTrain)-1; 

for i = 1:length(skipA) %can also change to length(input s times) = CHANGE BACK TO -1
Post_Stim_Matrix((1:100)+ newpairB(i),:) = outTrain((1:100)+ skipA(i),:); % 
end

trials = 1:60; 
times = 1:100; %100ms
%mean of all 
%make cells where each is one trial
rowDist = ones(1,ntrials)*100; % string of trials, each having 100 ms

%% make Post-stim times as a cell variable
Cell = mat2cell(Post_Stim_Matrix, rowDist); % each cell = one trial
electrodes = 1:60
%newly added - latency

latency_cell= Cell

for i=1:trials %NTRIALS

  l=latency_cell{i};
   
a(100,:) = 1;
latency_cell{i} = l
end


an = latency_cell{1};
bn = latency_cell{2};
cn = latency_cell{3};
dn = latency_cell{4};
en = latency_cell{5};
fn = latency_cell{6};
gn = latency_cell{7};
hn = latency_cell{8};
in = latency_cell{9};
jn = latency_cell{10};
Out = an+bn+cn +dn+en +fn +gn+hn+in+jn; % to get the first spike in first 10 trials 

First_spike = nan(1,electrode_no);

   for electrode_index = 1:60
 try
 First_spike(:, electrode_index) = find(Out(:, electrode_index), 1, 'first'); %fix
 catch
    First_spike = zeros(1,60)
 end
   end

   First_spike(First_spike ==100) = NaN;
  latencyHeatMap_MK(First_spike, channels,fileID)

   saveas(gcf,strcat(fileID,'_LatencyHM.png'));
try
   MEAN_fs = mean(First_spike)
catch
    MEAN_fs = 0
end

%
%
%% make a gif of over trials - some examples only
%{

for i = 1:length(C);

h =  figure;
    set(gcf,'Visible', 'off')
    imagesc(C{i});

    aesthetics
     axis tight manual 
% this ensures that getframe() returns a consistent size
filename = strcat(fileID, '.gif');
    ylabel('Time post stimulation (10 ms time bins)');
    xlabel('Electrode');
    cb = colorbar;
      
    % ylabel(cb, 'Spike count')
    ylabel(cb, 'Spike count');
    cb.TickDirection = 'out';
    set(gca,'TickDir','out');
    cb.Location = 'Southoutside';
    set(gca, 'FontSize', 14)
    xlimit_cbar = 4;
    caxis([0 4]);
    set(gcf,'Visible', 'off')
   drawnow

%capture as an image
    frame = getframe(h); 
    im = frame2im(frame); 
    [imind,cm] = rgb2ind(im,256); 


    if i == 1
          imwrite(imind,cm,filename,'gif', 'Loopcount',inf); 
      else 
          imwrite(imind,cm,filename,'gif','WriteMode','append') ;
          clear i 
    end 
end

%}


%% sum counts per trial
Post_stim_SCsum = cellfun(@sum,Cell,'uni',0); %sums SC for each trial
postStim_SC_Matrix = cell2mat(Post_stim_SCsum) ;%matrix of summed SC in 100 ms trials x electrodes 

ntrials = length(postStim_SC_Matrix); %GET RID OF -1
electrode_no = length(postStim_SC_Matrix');
%% Mean counts per electrode
Mean_SC = mean(postStim_SC_Matrix);
variance_SC = var(Mean_SC',0,1); %varaince between mean electrode spike counts
%within recordings
figure
bar(Mean_SC)
set(gcf,'Visible', 'off');
xlabel = 'Electrode';
ylabel('Mean Spike Count # 100 ms poststimulation');
ylim([0 25])
fileName = strcat(fileID, '_S_Counts_.mat');
saveas(gcf,strcat(fileName(1:end-4),'_Bargraph.png'));
MEAHeatMap_MK(Mean_SC,channels, 'counts');
saveas(gcf,strcat(fileID,'_stimSC_heatmap.png'));
SummeanNetwork_SC = sum(Mean_SC) ;   


a = sort(Mean_SC, 'descend');
top_10percent = a(1:6);
mean_SC_topNodes = mean(top_10percent);

%% Probablity of firing
BinaryM = logical(postStim_SC_Matrix); % binary matrix of wheather within one trial a spike occured 1x60
PPT_firing = sum(BinaryM)/ntrials ;%gives ppt of trials where spike was induced for each electrode
network_mean_P_firing = mean(PPT_firing)
top_10_prob = sort(PPT_firing, 'descend');
top_10_prob = top_10_prob(1:6);
mean_p_topNodes = mean(top_10_prob);

variance_PPT = var(PPT_firing', 0,1);
MEAHeatMap_MK(PPT_firing,channels, 'p')
saveas(gcf,strcat(fileID,'_FiringProbability_heatmap.png'));

responsive_network_size = sum(BinaryM); %number of electrodes 
network_ppt = sum(responsive_network_size>=1)/electrode_no; %ppt of the network which shows responses
%ms timewindow
save(strcat(fileID ,'_Response.mat') , 'PPT_firing', 'BinaryM', 'channels','Mean_SC','Post_Stim_Matrix', 'Cell',  'network_ppt', 'First_spike', 'SummeanNetwork_SC', 'variance_PPT', 'variance_SC', 'mean_SC_topNodes', 'network_mean_P_firing','mean_p_topNodes');
%change from baseline

end


%% PPT firing for each electrode paired graphs - OPTIONAL
files = dir('*OWT220208_1J*_cSpikes_L0_RP2*Response.mat')
 Matrix_PPT = zeros(5,60)
for i=1:length(files); %for the list of files
    fileID    =   files(i).name;  % at any one time the file in memory 
   %is selected from the list (i)
   data =  load(fileID);
   FR = data.PPT_firing
    Matrix_PPT(i, :) = FR(:,:)
end
 Matrix_PPT = Matrix_PPT'
  save(strcat(fileID(1:18),'_FiringRateE.mat'), 'Matrix_PPT' )

%repeat for spike sum

  Matrix_fr = zeros(5,60)
for i=1:length(files); %for the list of files
    fileID    =   files(i).name;  % at any one time the file in memory 
   %is selected from the list (i)
   data =  load(fileID);
   FR = data.Mean_SC
    Matrix_fr(i, :) = FR(:,:)
end
 Matrix_fr = Matrix_fr'
  save(strcat(fileID(1:18),'_SpikeCountsE.mat'), 'Matrix_fr' )
 %ORDER = BASELINE, HUB4, HUB6, PER4, PER6; unless some missing - just do
 %complete datasets for now

 %look at the effectof repeated trials

 %
files = dir('*OWT*RateE.mat')
for i = 1:length(files)
    fileID    =   files(i).name;  % at any one time the file in memory 
   %is selected from the list (i)
   data =  load(fileID);
   PPT = data.Matrix_PPT
   electrodes = [1:60]
electrodes= num2str(electrodes')
 labels = {'Baseline','Hub 4uA','Hub 6uA','Per 4uA', 'Per 6ua'};
h = figure
   parallelcoords(PPT, "Group", electrodes, Labels=labels)


legend('Location','northeastoutside', FontSize = 2)

saveas(gcf,strcat(fileID(1:18),'_PPTfiring_byelectrode.png'));
end
files = dir('*SpikeCountsE.mat')
for i = 1:length(files)
    fileID    =   files(i).name;  % at any one time the file in memory 
   %is selected from the list (i)
   data =  load(fileID);
   FR = data.Matrix_fr
   electrodes = [1:60]
electrodes= num2str(electrodes')
 labels = {'Baseline','Hub 4uA','Hub 6uA','Per 4uA', 'Per 6ua'};
h = figure
   parallelcoords(FR, "Group", electrodes, Labels=labels)


legend('Location','northeastoutside', FontSize=2)

saveas(gcf,strcat(fileID(1:18),'_SpikeCounts_byelectrode.png'));
end