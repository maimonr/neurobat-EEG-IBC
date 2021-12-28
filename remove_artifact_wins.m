function data = remove_artifact_wins(data,all_used_wins,winSize)
nWin = size(all_used_wins,1);
nChan = size(all_used_wins,2);

n_used_wins = size(data,1);

for ch_k = 1:nChan
    for w_k = 1:nWin
        if ~all_used_wins(w_k,ch_k)
            if nWin == n_used_wins
                data(w_k,ch_k,:,:) = NaN;
            else
                winIdx = 1 + (w_k-1)*winSize:(w_k*winSize);
                data(1,ch_k,:,winIdx) = NaN;
            end
        end
    end
end

end