function bpCorr = calculate_bandpower_corr(data,fs,freqBands,all_used_wins,winSize,varargin)

pnames = {'used_ch_pairs','winLength','overlapLength','bpMethod'};
dflts  = {'all',1,0.9,'windowed'};
[used_ch_pairs,winLength,overlapLength,bpMethod] = internal.stats.parseArgs(pnames,dflts,varargin{:});


nWin = size(data,1);
nChan = size(data,2)/2;
nT = size(data,3);
nFreq = size(freqBands,1);

switch bpMethod
    case 'windowed'
        
        data = reshape(data,[nWin nChan*2 1 nT]);
        data = remove_artifact_wins(data,all_used_wins,winSize);
        data = reshape(data,[nWin nChan*2 nT]);
        
        winLength = round(fs*winLength);
        winIdx = slidingWin(nT,winLength,round(fs*overlapLength));
        n_bp_win = size(winIdx,1);
        
        bandPower = nan(nWin,nChan*2,nFreq,n_bp_win);
        
    case 'hilbert'
        bandPower = nan(nWin,nChan*2,nFreq,nT);
end

for f_k = 1:nFreq
    [b,a] = butter(4,freqBands(f_k,:)/(fs/2),'bandpass');
    for win_k = 1:nWin
        for ch_k = 1:nChan*2
            chData = squeeze(data(win_k,ch_k,:));
            switch bpMethod
                case 'windowed'
                    bandPower(win_k,ch_k,f_k,:) = calculate_bp(chData,...
                        fs,freqBands(f_k,:),winIdx);
                case 'hilbert'
                    dataFilt = filtfilt(b,a,chData);
                    hilb = hilbert(dataFilt);
                    bandPower(win_k,ch_k,f_k,:) = abs(hilb)';
            end
        end
    end
end

if strcmp(bpMethod,'hilbert')
    bandPower = remove_artifact_wins(bandPower,all_used_wins,winSize);
end

switch used_ch_pairs
    case 'all'
        n_chan_pair = nChan.^2;
    case 'paired'
        n_chan_pair = nChan;
end

bpCorr = nan(nWin,nFreq,n_chan_pair);
for f_k = 1:nFreq
    for win_k = 1:nWin
        ch_k = 1;
        for ch1_k = 1:nChan
            for ch2_k = nChan+1:2*nChan
                if strcmp(used_ch_pairs,'paired') && (ch2_k ~= ch1_k + nChan)
                    continue
                end
                bpCorr(win_k,f_k,ch_k) = calculate_bpc(bandPower(win_k,[ch1_k,ch2_k],f_k,:));
                ch_k = ch_k + 1;
            end
        end
        
    end
end

end

function bpCorr = calculate_bpc(bandPower)

xBp = squeeze(bandPower(:,1,:,:));
yBp = squeeze(bandPower(:,2,:,:));
if ~all(isnan(xBp) | isnan(yBp),'all')
    nanIdx = ~any(isnan(xBp),2) & ~any(isnan(yBp),2);
    R = corrcoef(yBp(nanIdx),xBp(nanIdx));
    bpCorr = R(2);
else
    bpCorr = NaN;
end

end

function bp = calculate_bp(data,fs,freqBand,winIdx)

data = squeeze(data);
nWin = size(winIdx,1);

dataWin = data(winIdx);
usedIdx = ~any(isnan(dataWin),2);
bp = nan(1,nWin);
if any(usedIdx)
    bp(usedIdx) = bandpower(dataWin(usedIdx,:)',fs,freqBand);
end

end