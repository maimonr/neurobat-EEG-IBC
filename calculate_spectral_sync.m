function [ppc,pec,icoh,coh,wpli,pli] = calculate_spectral_sync(specData,varargin)

pnames = {'used_ch_pairs'};
dflts  = {'all'};
[used_ch_pairs] = internal.stats.parseArgs(pnames,dflts,varargin{:});

nFreq = size(specData,3);
nChan = size(specData,2)/2;
nWin = size(specData,1);

switch used_ch_pairs
    case 'all'
        n_chan_pair = nChan.^2;
    case 'paired'
        n_chan_pair = nChan;
end
[ppc,pec,icoh,coh,wpli,pli] = deal(nan(nWin,nFreq,n_chan_pair));
pcTransform = @(x,freq_k) log(abs(x(:,freq_k)).^2);
% pcTransform = @(x,freq_k) abs(x(:,freq_k));
for win_k = 1:nWin
    ch_k = 1;
    for ch1_k = 1:nChan
        for ch2_k = nChan+1:2*nChan
            if strcmp(used_ch_pairs,'paired') && (ch2_k ~= ch1_k + nChan)
                continue
            end
            if ~all(isnan(specData(win_k,ch1_k,:,:)) | isnan(specData(win_k,ch2_k,:,:)),'all')
                
                xCoef = squeeze(specData(win_k,ch1_k,:,:))';
                yCoef = squeeze(specData(win_k,ch2_k,:,:))';
                
                nanIdx = ~any(isnan(xCoef),2) & ~any(isnan(yCoef),2);
                xCoef = xCoef(nanIdx,:);
                yCoef = yCoef(nanIdx,:);
                
                [ppc(win_k,:,ch_k),pec(win_k,:,ch_k)] = calculate_power_corr(xCoef,yCoef,pcTransform);
                
                Sxy = xCoef.*conj(yCoef);
                [coh(win_k,:,ch_k),icoh(win_k,:,ch_k)] = calculate_coh(xCoef,yCoef,Sxy);
                
                [pli(win_k,:,ch_k),wpli(win_k,:,ch_k)] = calculate_wpli(Sxy);
                
            end
            ch_k = ch_k + 1;
        end
    end
end
end

function [ppc,pec] = calculate_power_corr(xCoef,yCoef,pcTransform)

Yprojx = imag(yCoef.*(conj(xCoef)./abs(xCoef)));
Xprojy = imag(xCoef.*(conj(yCoef)./abs(yCoef)));

nFreq = size(xCoef,2);
Rppc = nan(nFreq,2);
Rpec = nan(nFreq,1);

for freq_k = 1:nFreq
    Rppc(freq_k,1) = corr(pcTransform(Yprojx,freq_k),pcTransform(xCoef,freq_k));
    Rppc(freq_k,2) = corr(pcTransform(Xprojy,freq_k),pcTransform(yCoef,freq_k));
    
    Rpec(freq_k) = corr(pcTransform(yCoef,freq_k),pcTransform(xCoef,freq_k));
end

ppc = mean(Rppc,2);
pec = Rpec;


end

function [coh,icoh] = calculate_coh(xCoef,yCoef,Sxy)

Sxy = nansum(Sxy);
Sxx = nansum(xCoef.*conj(xCoef));
Syy = nansum(yCoef.*conj(yCoef));

Cxy = Sxy./(sqrt(Sxx).*sqrt(Syy));

coh = abs(Cxy);
icoh = abs(imag(Cxy));

end

function [pli,wpli] = calculate_wpli(Sxy)

Ixy = imag(Sxy);
pli = abs(mean(sign(Ixy)));
wpli = (abs(mean(Ixy)))./mean(abs(Ixy));

end