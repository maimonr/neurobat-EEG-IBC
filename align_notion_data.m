function [dataAligned,tAligned] = align_notion_data(t,data)

fs = 250;
maxDT = (1/fs)*4;

nDevice = size(data,2);
nFile = size(data,1);

tAligned = cell(nFile,1);
dataAligned = cell(nFile,nDevice);

for f_k = 1:nFile
    
    tCurrent = t(f_k,:);
    dataCurrent = data(f_k,:);
    discontIdx = cellfun(@(t) find(abs(diff(t)) > maxDT,1,'first'),tCurrent,'un',0);
    
    for d_k = 1:nDevice
        if ~isempty(discontIdx{d_k})
            assert(discontIdx{d_k}/length(tCurrent{d_k}) > 0.95)
            tCurrent{d_k} = tCurrent{d_k}(1:discontIdx{d_k});
            dataCurrent{d_k} = dataCurrent{d_k}(1:discontIdx{d_k},:);
        end
    end
    
    assert(all(cellfun(@(t) all(diff(t) <   maxDT & diff(t)>= -1/fs),tCurrent)))
    
    for d_k = 1:nDevice
       x = [1:length(tCurrent{d_k})]';
       p = polyfit(x,tCurrent{d_k},1);
       y = polyval(p,x);
       assert(all(abs(y - tCurrent{d_k}) < maxDT/2))
       tCurrent{d_k} = y;
    end
    
    tRange = [max(cellfun(@min,tCurrent)) min(cellfun(@max,tCurrent))];
    assert(all(cellfun(@(t) sum(t<tRange(1) | t>tRange(2))/length(t) < 0.95,tCurrent)))
    
    dataCurrent = cellfun(@(d,t) d(t>=tRange(1) & t<=tRange(2),:),dataCurrent,tCurrent,'un',0);
    tCurrent = cellfun(@(t) t(t>=tRange(1) & t<=tRange(2)),tCurrent,'un',0);
    
    tAligned{f_k} = tCurrent{1};
    
    nanIdx = cellfun(@(d) false(size(d)),dataCurrent,'un',0);
    for d_k = 1:nDevice
        if d_k == 1
            dataCurrent{d_k} = dataCurrent{d_k};
        else
            dataCurrent{d_k} = interp1(tCurrent{d_k},dataCurrent{d_k},tAligned{f_k});
            nanIdx{d_k} = isnan(dataCurrent{d_k});
            assert(sum(any(nanIdx{d_k},2))/size(nanIdx{d_k},1) < 0.95)
        end
    end
    
    nanIdx = any(cat(3,nanIdx{:}),[2 3]);
    
    dataAligned(f_k,:) = cellfun(@(d) d(~nanIdx,:),dataCurrent,'un',0);
    tAligned{f_k} = tAligned{f_k}(~nanIdx);
    
end

end