function results = run_notion_eeg_corr(expDate)
expTypes = {'gaze','back','close'};
% expTypes = {'gaze','mirror','back','close'};

selfCorr = false;
self_corr_device = 1;

keepWins = false;

freqBands = [1 4; 4 8; 8 12; 13 30];

fs = 250;
winLength = 2*fs;
spec_win_length = 1;
overlap = 0;

nExp = length(expTypes);
results = cell(1,nExp);
tt = tic;
for exp_k = 1:nExp
    
    [data, t] = read_notion_json(expTypes{exp_k},expDate);
    
    if selfCorr
        data(:,setdiff(1:2,self_corr_device)) = data(:,self_corr_device);
        t(:,setdiff(1:2,self_corr_device)) = t(:,self_corr_device);
    end
    
    dataAligned = align_notion_data(t,data);
    [dataWin,usedWins] = window_eeg_data(dataAligned,winLength,overlap);
    
    [allData, all_used_wins] = collapse_windowed_data(dataWin,usedWins,keepWins);
    
    bpCorr = calculate_bandpower_corr(allData,fs,freqBands,all_used_wins,winLength);
    
    specData = calculate_eeg_spec(allData,spec_win_length);
    specData = remove_artifact_wins(specData,all_used_wins,winLength);
    
    [ppc,pec] = calculate_ppc(specData);
    [icoh,coh] = calculate_icoh(specData);
    
    currentResults = struct;
    currentResults.ppc = ppc;
    currentResults.pec = pec;
    currentResults.icoh = icoh;
    currentResults.coh = coh;
    currentResults.bp = bpCorr;
    results{exp_k} = currentResults;
    
    fprintf('%s finished, %f s elapsed\n',expTypes{exp_k},toc(tt));

end

results = cell2struct(results,expTypes,2);

end