function [allData, all_used_wins] = collapse_windowed_data(dataWin,usedWins,varargin)

pnames = {'keepWins','replaceArtifacts'};
dflts  = {false,false};
[keepWins,replaceArtifacts] = internal.stats.parseArgs(pnames,dflts,varargin{:});


nFile = size(dataWin,1);
[all_used_wins,allData] = deal(cell(1,nFile));
for f_k = 1:nFile
    allData{f_k} = cat(3,dataWin{f_k,:});
    allData{f_k} = num2cell(allData{f_k},[1 3]);
    allData{f_k} = cellfun(@(y) squeeze(y)',allData{f_k},'un',0);
    all_used_wins{f_k} = cat(2,usedWins{f_k,:});
end
all_used_wins = vertcat(all_used_wins{:});

allData = [allData{:}];

if replaceArtifacts
    for k = 1:length(allData)
        allData{k}(~all_used_wins(k,:),:) = NaN;
    end
end

if keepWins
    allData = cat(3,allData{:});
else
    allData = [allData{:}];
end

allData = permute(allData,[3 1 2]);

end
