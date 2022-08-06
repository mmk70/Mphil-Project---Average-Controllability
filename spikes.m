%% set up parameters -change accordingly
clear all
pwd()
HomeDir = 'C:\Users\user\OneDrive - University of Cambridge\Documents\MATLAB\Project'; %set here source data
RecordingDir ='D:'
cd(RecordingDir)
addpath(genpath(HomeDir))
fs = 25000; %sampling rate /s
sampling_rate= 25000 %sampling rate /s (some functions require sampling_rate instead of fs)
duration = round(300 * fs); %adjust time if needed, currently set to 300 s
start_time = 0 ;
end_time = 300; %s
num_channels = 60 ;
cd(RecordingDir)


%% get spikes - fast method (option three)

files =  dir('MEC*RP2.mat')
for i=1:length(files); %for the list of files
    fileID    =   files(i).name;  % at any one time the file in memory 
   %is selected from the list (i)
   data =  load(fileID);
 
%SPIKE DETECTION METHODS AND PARAMETERS
%option one: two methods, three parameters for each
%meths   ={'Manuel';'cwt'};
%params  =[4,5,6;                                % Threshold
 %   -0.1254,0,0.1254]';                         % L parameter

%option two: two methods, one parameter for each
% meths   ={'Manuel';'cwt'};
% params  =[5;0]';                               % threshold; L parameter

%option three: one method, one parameter
meths   ={'cwt'};                              % one method only for speed
params  =[4,0]; %threshold, L parameter - CHANGED FROM 0 TO 0.12 and back

%option four: one method, three parameters
% meths   ={'Manuel'};
% params  =[4,4.25,4.5;                                % Threshold
%     -0.1254,0,0.1254]';                         % L parameter
%


%% get spikes

refPeriod_ms = 2.0; %choose refractory period in ms 
for m = 1:length(meths)
    method              =       meths{m}
    
    for p = 1:size(params,1);
        L               =       params(p+size(params,1))
        multiplier      =       params(p)
        
        batchGetSpike_fcn(RecordingDir,files,method,multiplier,L,refPeriod_ms) %may have to fiddle with that
        progressbar([],p/size(params,1),[]);       %update parameter progress
    end
   progressbar([],[],m/length(meths));            %update % methods done
end
warning('on','MATLAB:load:variableNotFound')
end






%% get Spike Matrix
files = dir('*OWT*Spikes*RP2.mat');
fs =25000
for i=1:length(files);
    fileName    =   files(i).name;
    
    data = load(fileName)
channels = data.channels;
spikeMatrix = full(data.cSpikes);

   %% downsampling BY 100 MS
   num_channels = 60
     recordDuration = round(length(spikeMatrix)); %in samples
    % recordDuration_100ms = 3000
    downSpikeMatrix = downSampleSum_org(spikeMatrix, recordDuration * 1/fs);
    downFactor = 25000; %same as fs
     downSpikeMatrix = downSampleSum(spikeMatrix, recordDuration/downFactor);
    %downFactor_100ms = 2500
   downTrain = reshape(spikeMatrix, [], recordDuration, num_channels);
    downTrain = sum(downTrain, 1);
    outTrain = reshape(downTrain, recordDuration, num_channels);

    new_fs = fs/downFactor;

%% raster plot
 figure
   imagesc(log10(outTrain'))
  imagesc(downSpikeMatrix')
    
  aesthetics
 ylabel('Electrode')
 xlabel('Time (s)')
 cb = colorbar;
 ylabel(cb, 'Spike count')
 ylabel(cb, 'Firing Rate (Hz)')
   cb.TickDirection = 'out';
 cb.Ticks = 0:5; % for slice 5 specifically
 set(gca,'TickDir','out');
 cb.Location = 'Southoutside';
  cb.Box = 'off';
 set(gca, 'FontSize', 14)
 ylimit_cbar = 3;
 caxis([0,ylimit_cbar]) %set colorbar axis limits; also adjusts colour
below command does not adjust colour hence need for caxis command
remove decimal ticks e.g. 0.5 Hz
   cb.Ticks = linspace(0,ylimit_cbar,ylimit_cbar/1+1);%(start,end,number of numbers)
  %below is for plotting log scale, labelling in raw spike count
  cb.TickLabels = 10.^(linspace(0,ylimit_cbar,ylimit_cbar+1));
    
  yticks([1, 10:10:60])
    
  fileName1=files(i).name;
if  contains(files(i).name,'_') %remove underscores for title
  fileName1(strfind(fileName1,'_'))=' ';
fileName1=strcat('{',fileName1,'}');
 title({strcat(fileName1(1:end-4),' Raster'),' '});
else
   title({strcat(files(i).name(1:end-4),' Raster'),' '});
end
  %make title font smaller
  ax = gca;
 ax.TitleFontSizeMultiplier = 0.7;
    
save raster as PNG
fprintf(strcat('\n','\n',files(i).name(1:end-4),' saving raster...', '\n','\n'))
saveas(gcf,strcat(fileName(1:end-4),'_Raster.png'));
 


%% get firing rates
SpikeCounts = sum(spikeMatrix);
FiringRates=SpikeCounts/(length(spikeMatrix)/fs);


%% get burst info
samplingRate=fs;
method ='Bakkum';
% Set N = 30 (min number of bursts)
%ensure bursts are excluded if fewer than 3 channels (see inside burst Detect
%function)
%to change min channels change line 207 of burstDetect.m
%to change N (min n spikes) see line 170 of burstDetect.m
N = 30; %min number of bursts (set to 30)
minChannel = 3;%no of channels participating, 3
[burstMatrix, burstTimes, burstChannels] = burstDetect(spikeMatrix, method, samplingRate,N,minChannel);
nBursts=length(burstMatrix);
BurstRate = nBursts/(length(spikeMatrix(:,1))/fs)
BurstchNo =  zeros(1,nBursts);
for i = 1:length(burstChannels)
BurstchNo(i) = max(size(burstChannels{i}));
end
Output.channels = channels
 Output.maxnetworksize = max(BurstchNo);
 Output.averageChNo = mean(BurstchNo);
 Output.meannetwroksize = min(BurstchNo);
 Output.sdNs = std(BurstchNo);
 Output.BurstRate = BurstRate;
 Output.SpikeCounts = SpikeCounts;
Output.FiringRates = FiringRates;
save(strcat(fileName,'_BurstOutput.mat'), "Output" )

end
%% make MEA heatmaps
option = 'rate';
batch_getHeatMaps_fcn(files,option)


%% get Adj M
files = dir('OWT*cSpikes*RP2.mat')
parallel = 1;
delta_t = 0.01; % in s
thr = 0.1;
desired_hz = 1000;
duration_s =300
method = 'tileCoef'
sync_win_s = 0.05; %synchroncity window in seconds; e.g. 1 is +/- 1 s is considered synchronous (2DeltaT = 2 s)
rec_length_s = 300;
fs = 25000;
rec_length_samps = fs * rec_length_s;

%% downsampling:
num_samples = 1000 * rec_length_s; %defaults to 1 sample per ms
ds_rate = num_samples/rec_length_samps; %scalar for time_window %downsampling rate
sync_win = sync_win_s * ds_rate; %downsample the sync_win

batch_getAdj_fcn(method,files,sync_win,num_samples,ds_rate);
%add image of adjM


%% AC node detection - load individually 
data= load('MEC220425_6D_DIV53_BASELINE_cSpikes_L0_RP2_adjM_0.05.mat')
 %load AdjM for Baseline recordings
 channels = data.channels;
 adjM = data.adjM;
 AC = ave_control(adjM) %calculate average controllablity for each node
[A,I] = sort(AC,'descend') %sort nodes by AC
channels_sorted = channels(I); %sort corresponding channel number by AC
channels_sorted(1) %get top AC node to stimulate
A(1) %get value of top AC node

X= load('MEC220425_6D_DIV53_BASELINE_cSpikes_L0_RP2.mat_BurstOutput.mat') %load firing rates for Baseline recording
x = X.Output
x.FiringRates(I) %find INDEX lowest AC active node (0.015>) 
channels_sorted(53) %get corresponding channel number FROM LINE ABOVE
A(53) %get AC value oflowest AC node
