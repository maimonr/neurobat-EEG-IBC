function [dataWin,usedWins] = window_eeg_data(data,winLength,varargin)

pnames = {'artifact_quantiles','artifactThresh'};
dflts  = {[0.995 0.005],0.95};
[artifact_quantiles,artifactThresh] = internal.stats.parseArgs(pnames,dflts,varargin{:});

nFile = size(data,1);
nDevice = size(data,2);

dataWin = get_windowed_eeg(data,winLength);

usedWins = deal(cell(1,nDevice));
for d_k = 1:nDevice
    allData = cat(1,data{:,d_k});
    [all_data_z,mu,sigma] = zscore(allData);
    mu = reshape(mu,[1 1 length(mu)]);
    sigma = reshape(sigma,[1 1 length(sigma)]);
    q = quantile(all_data_z(:),artifact_quantiles);
    %     q = nan(2,1,nChan);
    %     q(:,1,:) = quantile(allData,artifact_quantiles,1);
    
    for f_k = 1:nFile
        Z = (dataWin{f_k,d_k} - mu)./sigma;
        artifactWin = Z < q(1) & Z > q(2);
        usedWins{f_k,d_k} = squeeze(sum(artifactWin,1))/size(artifactWin,1) > artifactThresh;
    end
end


end

function dataWin = get_windowed_eeg(data,winLength)

nFile = size(data,1);
nDevice = size(data,2);
nChan = size(data{1},2);

winIdx = cellfun(@(d) slidingWin(size(d,1),winLength,0),data(:,1),'un',0);
dataWin = cell(size(data));
for d_k = 1:nDevice
    for f_k = 1:nFile
        dataWin{f_k,d_k} = nan(size(winIdx{f_k},1),winLength,nChan);
        for ch_k = 1:nChan
            dataTmp = data{f_k,d_k}(:,ch_k);
            dataWin{f_k,d_k}(:,:,ch_k) = dataTmp(winIdx{f_k});
        end
    end
end
dataWin = cellfun(@(d) permute(d,[2 1 3]),dataWin,'un',0);

end