%% metrics - get metrics from AdjM analysis of spontaneous recordings
%% Metrics output

%AM = output file
files = dir('*adjM_0.05.mat')

for i = 1:length(files)
    fileID= files(i).name
    data = load(fileID)
    adjM = data.adjM
adjM(isnan(adjM))=0
channels = data.channels
   AM.channels = data.channels

%average controllability
data = load('OWT220207_2D_DIV63_BASELINE_cSpikes_L0_RP2_adjM_0.05.mat')
adjM = data.adjM
channels = data.channels
AC = ave_control(adjM);
MEAHeatMap_MK(AC,channels, 'AC'); %get heatmap

% modal controllablity
AM.MC = modal_control(adjM,'complex');

% participation coef
M     = community_louvain(adjM);
 AM.P = participation_coef(adjM, M);

% efficiency 
  AM.LE = efficiency_wei(adjM, 2);

% betweenness c
AM.BC = betweenness_wei(adjM);

% node strength
AM.NS = sum(adjM,2); %sum of adj m row  (row 1 = node 1 containing sttc ie edgeweights)
AM.meanNS = mean(AM.NS)
AM.meanAC = mean(AM.AC)
AM.meanBC = mean(AM.BC)
AM.meanP = mean(AM.P)
AM.meanLE = mean(AM.LE)

AM.zscoreAC = zscore(AM.AC)
AM.zscoreMC = zscore(AM.MC)
AM.zscoreBC = zscore(AM.BC)
AM.zscoreLE = zscore(AM.LE)
AM.zscorePC = zscore(AM.P)
AM.zscoreNS = zscore(AM.NS)
Hub_matrix = [AM.zscoreNS AM.zscoreBC AM.zscoreLE AM.zscorePC];
AM.meanHubScore = mean(Hub_matrix,2); %mean of all
AM.adjM = adjM

[A,ACIndex] = sort(AM.AC,'descend');
AM.ACchannels_sorted = channels(ACIndex)
[H,HIndex] = sort(AM.meanHubScore,'descend') 
AM.HUBchannels_sorted = channels(HIndex)
[N,NIndex] = sort(AM.NS, 'descend') 
AM.NSchannels_sorted = channels(NIndex)
[N,MCIndex] = sort(AM.MC, 'descend') 
AM.MCchannels_sorted = channels(MCIndex)
save(strcat(fileID,'_RecordingMetrics.mat'), 'AM' )
end




