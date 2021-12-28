function results = run_eeg_corr(expDate,varargin)

pnames = {'used_ch_pairs','offset_to_drop','keepWins','specMethod',...
    'winLength','fs','specParams','freqLims','freqBands'};
dflts  = {'all',5,false,'cwt',...
    5,128,[],[1 40],[]};
[used_ch_pairs,offset_to_drop,keepWins,specMethod,...
    winLength,fs,specParams,freqLims,freqBands] = internal.stats.parseArgs(pnames,dflts,varargin{:});

switch datestr(expDate,'yyyymmdd')
    case '20211210'
        expTypes = {'mutual_gaze','eyes_closed','visual_flicker_20hz','finger_tapping',...
            'metronome_180bpm','finger_tapping_metronome_180bpm'};
    case '20211206'
        expTypes = {'mutual_gaze','eyes_closed','drifting_gratings','finger_tapping'};
end

if isempty(freqBands)
    freqBands = [1 4; 4 8; 8 12; 13 30];
end

winLength = winLength*fs;

if isempty(specParams)
    specParams = struct('spec_win_length',winLength,'timeBandwidth',60,...
        'freqLims',freqLims);
end
specParams.specMethod = specMethod;

nExp = length(expTypes);
results = cell(1,nExp);
tt = tic;
for exp_k = 1:nExp
    
    [data, t, deviceNames, chNames] = read_emotiv_data(expDate,expTypes{exp_k},...
        'offset_to_drop',offset_to_drop);
    ch_pair_labels = get_ch_pair_labels(deviceNames,chNames,used_ch_pairs);
    
    dataAligned = align_emotiv_data(t,data);
    [dataWin,usedWins] = window_eeg_data(dataAligned,winLength);
    
    [allData, all_used_wins] = collapse_windowed_data(dataWin,usedWins,'keepWins',keepWins);
    
    [specData,specFreqs] = calculate_eeg_spec(allData,fs,specParams);
    specData = remove_artifact_wins(specData,all_used_wins,winLength);
    
    currentResults = struct;
    
    [ppc,pec,icoh,coh,wpli,pli] = calculate_spectral_sync(specData,'used_ch_pairs',used_ch_pairs);
    bpCorr = calculate_bandpower_corr(allData,fs,freqBands,all_used_wins,winLength,...
        'used_ch_pairs',used_ch_pairs);
    
    currentResults.ppc = ppc;
    currentResults.pec = pec;
    currentResults.icoh = icoh;
    currentResults.coh = coh;
    currentResults.pli = pli;
    currentResults.wpli = wpli;
    currentResults.bp = bpCorr;
    currentResults.specFreqs = specFreqs;
    currentResults.pairNames = {ch_pair_labels};
    
    results{exp_k} = currentResults;
    
    fprintf('%s finished, %f s elapsed\n',expTypes{exp_k},toc(tt));
    
end

results = cell2struct(results,expTypes,2);
paramValues = {used_ch_pairs,offset_to_drop,keepWins,specMethod,...
    winLength,fs,specParams,freqLims,freqBands};
results.params = cell2struct(paramValues,pnames,2);

end

function ch_pair_labels = get_ch_pair_labels(deviceNames,chNames,used_ch_pairs)

chLabels = cell(1,length(deviceNames));
for d_k = 1:length(deviceNames)
    chLabels{d_k} = cellfun(@(ch) [deviceNames{d_k} '_ch_' ch],chNames,'un',0);
end

nChan = length(chNames);
switch used_ch_pairs
    case 'all'
        n_chan_pair = nChan.^2;
    case 'paired'
        n_chan_pair = nChan;
end

ch_pair_labels = cell(n_chan_pair,1);

ch_k = 1;

for ch1_k = 1:nChan
    for ch2_k = 1:nChan
        if strcmp(used_ch_pairs,'paired') && ch1_k ~= ch2_k
            continue 
        end
        ch_pair_labels{ch_k} = strjoin({chLabels{1}{ch1_k},chLabels{2}{ch2_k}},'-');
        ch_k = ch_k + 1;
    end
end

end