%Remove artefact
%identify artefact time BY SYNCHRONY OF OCCURANCE at EARLIEST TIMEPOINT
files = dir('*cSpikes*RP2.mat');

for i = 1:length(files)
    fileID= files(i).name;
    data = load(fileID);
 
    spikes = data.cSpikes;
    channels = data.channels;
    spikes = full(spikes);
    
    electrode_no = 60;
  z = reshape(spikes, [], length(spikes)/2, electrode_no);
    zz = sum(z, 1); 
  Condensed_spikes = reshape(zz, length(spikes)/2, electrode_no);
sumOfSpikes_PerRow = sum(Condensed_spikes,2); %GET Sum of spikes per row to find 
%where the most synchrony occurs = artefact
Trial_no = 60 %change accordingly
 [M,I] = maxk(sumOfSpikes_PerRow,Trial_no); % find index of highest synchrony x 60 (no of trials)
 if M(1) >40 %SET TO 40 NODES
I = I *2; % artefact time in frame
Indexdelay = I +1; %some happen one frame after
Indexpre = I -1;
Indexpre2 = I -2;
Indexpre3 = I -3;
%zero out 4.5ms before artefact for sure arrives
 Artefact_times = [I; Indexdelay; Indexpre; Indexpre2; Indexpre3]; %in frames
 Artefact_times(Artefact_times<=0) = []; %get rid of zeros
 
spikes(Artefact_times,:) = 0;
spikes = sparse(spikes);
%save as data with no artefacts
save(strcat(fileID, '_Nortefact.mat'), 'spikes', 'channels');
 else
     disp('No artefact detected')
 end
end

%LATENCY: do the first few trials, maybe even just the first one idk 
%should we just do the time of first spike in filtered data?

