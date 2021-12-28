function plot_ibc_results(res)
freqBands = [1 4; 4 8; 8 12; 13 30];
freqNames = {'delta','theta','alpha','beta'};
expTypes = {'mutual_gaze','eyes_closed','visual_flicker_20hz','finger_tapping',...
    'metronome_180bpm','finger_tapping_metronome_180bpm'};
expStrs = {'mg','ec','vf','ft','met','met_ft'};


avgType = 'window';

nExp = length(expTypes);
nF = size(freqBands,1);
clf
usedMeasures = {'ppc','pec','icoh','bp'};
tiledlayout(length(usedMeasures),nF,'TileSpacing','compact');

for meas = usedMeasures
    nChan = size(res.(expTypes{1}).(meas{1}),3);
    nWin = max(cellfun(@(exp) size(res.(exp).(meas{1}),1),expTypes));
    
    switch avgType
        case 'channel'
            allMeas = nan(nWin,nF,nExp);
            avgIdx = 3;
            usedIdx = 1;
        case 'window'
            allMeas = nan(nChan,nF,nExp);
            avgIdx = 1;
            usedIdx = 3;
    end
    
    for freq_k = 1:size(freqBands,1)
        for exp_k = 1:nExp
            n_exp_samples = size(res.(expTypes{exp_k}).(meas{1}),usedIdx);
            if ~strcmp(meas{1},'bp')
                idx = freqBands(freq_k,1):freqBands(freq_k,2);
                allMeas(1:n_exp_samples,freq_k,exp_k) = ...
                    squeeze(mean(res.(expTypes{exp_k}).(meas{1})(:,idx,:),[2 avgIdx],'omitnan'));
            else
                allMeas(1:n_exp_samples,freq_k,exp_k) = ...
                    squeeze(mean(res.(expTypes{exp_k}).(meas{1})(:,freq_k,:),avgIdx,'omitnan'));
            end
        end
    end

    for freq_k = 1:nF
        nexttile
        boxplot(squeeze(allMeas(:,freq_k,:)))
        h = gca;
        h.XTickLabels = expStrs;
        title([meas{1} ' - ' freqNames{freq_k}])
        axis square
    end
end
end