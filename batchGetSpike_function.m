function batchGetSpike_function(datadir,files,method,multiplier,L,refPeriod_ms,my_email)
% read all .mat files in directory, look for spikeMatrix or dat
% then extract spikes, save as SPARSE MATRIX


% note that this overwrites rather than append to the info file

% I initially wanted to do this together with batch processing,
% but since I am still experimenting on the results that are obtained
% from the various spike detection algorithms, I want to save them as
% sparse spike file first, so that in the future I can play around with
% them without having to load the raw data again (unless I want to try
% out another spike detection parameter / method)



%% some parameters

%files = dir('190830_slice1stim5.mat');  % where your .mat files are
%files = files(~contains({files.name}, 'Spikes'));
%files = files(~contains({files.name}, 'Filt'));%remove unwanted files
%disp(strcat('number of files to do = ',num2str(length(files)))
% variable name containing MEA voltage recording should be either of
% these two:
voltageVar2 = 'dat';
voltageVar = 'electrodeMatrix';
% assume it takes the form numSamp x numChannels
% samplingRate = 25000;
%progressbar
for file = 1:length(files)
    if strcmp(method,'cwt')
        fileName = strcat(files(file).name(1:end-4), '_cSpikes_L',num2str(L),'_RP',num2str(refPeriod_ms), '.mat');
    elseif strcmp(method,'Manuel')
        fileName = strcat(files(file).name(1:end-4), '_mSpikes_',num2str(multiplier),'_RP',num2str(refPeriod_ms), '.mat');
    elseif strcmp(method,'abs')
        fileName = strcat(files(file).name(1:end-4), '_aSpikes_',num2str(multiplier),'_RP',num2str(refPeriod_ms), '.mat');
    end
    if ~exist(fileName)
%         try
            try
                data = load(files(file).name, voltageVar);
                data = data.(voltageVar);
                channels = load(files(file).name, 'channels');
                channels = channels.('channels');
                fs = load(files(file).name, 'fs');
                fs = fs.('fs');
            catch
                data = load(files(file).name, voltageVar2);
                data = data.(voltageVar2);
                channels = load(files(file).name, 'channels');
                channels = channels.('channels');
                fs = load(files(file).name, 'fs');
                fs = fs.('fs');
                fprintf('Data loaded successfully \n')
            end
             %data = data.(voltageVar); % since matlab load struct
             %data = electrodeMatrix
            %% detect spikes
            % tic;
            % mSpikes = sparse(getSpikeMatrix(data, 'Manuel', 5));
            % tSpikes = sparse(getSpikeMatrixAlex(data, 'Tim', 14)); %AD edited to use my getSpikeMatrix.m edited script; threshold changed from 8 to 12
            % pSpikes = sparse(getSpikeMatrix(data, 'Prez', 4));
            %L=-0.1254; %changed L to 0 to confirm it changes spike rates - it did indeed increase them as it should
            % L = log(false detection / false omission) / 36.7368
            % e.g. if you want low sensitivity (cost of false detection 1000x
            % greater), do log(1000)/36.7368 to derive L
            
            %% get spikes and save
   
            if strcmp(method,'cwt')
                fileName = strcat(files(file).name(1:end-4), '_cSpikes_L',num2str(L),'_RP',num2str(refPeriod_ms), '.mat');
                if ~exist(fileName)
                    [spikeMatrix,~,thresholds] = getSpikeMatrixAlex_fcn(data, 'cwt', multiplier, L,fs,refPeriod_ms); %AD added, need to adjust loss parameter
                    spikeMatrix(:,find(channels == 15)) = zeros(size(spikeMatrix(:,find(channels == 15)))); %remove spikes from ref channel
                    cSpikes = sparse(spikeMatrix);
                    %toc
                    % save(fileName, 'mSpikes', 'tSpikes', 'pSpikes');
                    save(fileName, 'cSpikes','channels','thresholds','refPeriod_ms','method','fs','multiplier');
                else
                    disp('spikes already detected')
             
%        % catch
%             %error_email_fcn(my_email,files(file).name);
%        % end
    
        %         disp('spikes already detected')
    end
    progressbar(file/length(files),[],[]);
    %D:\MECP2_2019_AD\Scripts_and_Output\S2.1.SpikeMatrix_Scripts
    %this is where current script is located and called functions are
            end
    end
end
end


