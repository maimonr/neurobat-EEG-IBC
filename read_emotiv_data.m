function [data, t, deviceNames, chNames] = read_emotiv_data(expDate,expType,varargin)

pnames = {'fs','filterFreqs','offset_to_drop'};
dflts  = {128,[0.5 50],5};
[fs,filterFreqs,offset_to_drop] = internal.stats.parseArgs(pnames,dflts,varargin{:});

[b,a] = butter(7,filterFreqs/(fs/2),'bandpass');
exp_day_str = datestr(expDate,'yyyymmdd');

baseDir = 'C:\Users\BatLab\Documents\eeg_data\';
eegData = load(fullfile(baseDir,['emotiv_data_' exp_day_str '.mat']),expType,'ids');
deviceNames = eegData.ids;
eegData = eegData.(expType);

used_ch_names = {'AF3';'F7';'F3';'FC5';'T7';'P7';'O1';'O2';'P8';'T8';'FC6';'F4';'F8';'AF4'};
all_ch_names = cellstr(squeeze(eegData{1}.chNames));
chIdx = ismember(all_ch_names,used_ch_names);
chNames = all_ch_names(chIdx);

deviceNames = cellstr(deviceNames);
[nFile,nDevice] = size(eegData);
[data,t] = deal(cell(nFile,nDevice));

for d_k = 1:nDevice
    
    for f_k = 1:nFile
        data{f_k,d_k} = filtfilt(b,a,eegData{f_k,d_k}.data(:,chIdx));
        t{f_k,d_k} = eegData{f_k,d_k}.timestamps;
    end
    
end

if offset_to_drop > 0
    
    n_sample_to_drop = offset_to_drop*fs;
    data = cellfun(@(d) d(n_sample_to_drop:end-n_sample_to_drop,:),data,'un',0);
    t = cellfun(@(d) d(n_sample_to_drop:end-n_sample_to_drop),t,'un',0);
    
end

end