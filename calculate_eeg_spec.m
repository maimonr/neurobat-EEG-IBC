function [specData,f] = calculate_eeg_spec(allData, fs, specParams)

switch specParams.specMethod
    
    case 'ft_mtmconvol'
        
        [specData,f] = calculate_ft_eeg_spec(allData,fs,specParams);
        
    case 'cwt'
        
        [specData,f] = calculate_cwt_eeg_spec(allData,fs,specParams);
        
        
end

end

function [specData,f] = calculate_cwt_eeg_spec(allData,fs,specParams)

nWin = size(allData,1);
nChan = size(allData,2);
nT = size(allData,3);

fb = cwtfilterbank('SamplingFrequency',fs,'FrequencyLimits',specParams.freqLims,...
    'SignalLength',nT,'TimeBandwidth',specParams.timeBandwidth);

f = fb.centerFrequencies;
nF = length(f);

specData = nan(nWin,nChan,nF,nT);

for win_k = 1:nWin
   for ch_k = 1:nChan
       specData(win_k,ch_k,:,:) = cwt(squeeze(allData(win_k,ch_k,:)),'FilterBank',fb);
   end
end

end

function [specData,f] = calculate_ft_eeg_spec(allData,fs,specParams)

cfg              = [];
cfg.output       = 'fourier';
cfg.channel      = 'all';

cfg.pad          = 'nextpow2';
cfg.method       = 'mtmconvol';
cfg.taper        = 'hanning';
cfg.foi          = specParams.freqLims(1):1:specParams.freqLims(2);
cfg.t_ftimwin    = ones(length(cfg.foi),1).*specParams.specWin;   % length of time window = specWin
cfg.toi          = 'all';

dT = 1/fs;

t = 0:dT:(size(allData,3)/fs)-dT;

labels = cellfun(@(ch_k) ['ch_' num2str(ch_k)],num2cell(1:size(allData,2)),'un',0);

nWin = size(allData,1);

ft_T = repmat({t},1,nWin);
ft_data = num2cell(allData,[2 3]);
ft_data = cellfun(@squeeze,ft_data,'un',0);

ftData = struct('label',{labels'},'fsample',fs,'trial',{ft_data},'time',{ft_T});
tfOut = ft_freqanalysis(cfg,ftData);

specData = tfOut.fourierspctrm;
f = cfg.foi;

end

