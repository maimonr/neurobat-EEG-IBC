function plot_spectral_ibc_results(results,usedMeasure)
expTypes = {'mutual_gaze','eyes_closed','visual_flicker_20hz','finger_tapping',...
    'metronome_180bpm','finger_tapping_metronome_180bpm'};

use_matched_pairs = false;
pairNames = results.(expTypes{1}).pairNames{1};
chNames = cellfun(@(pName) cellfun(@(x) x(3), cellfun(@(s) strsplit(s,'_'), strsplit(pName,'-'),'un',0)),pairNames,'un',0);
matched_pair_idx = cellfun(@(ch) strcmp(ch{:}),chNames);
if use_matched_pairs
    chIdx = matched_pair_idx;
else
    chIdx = true(size(matched_pair_idx));
end

figure
tiledlayout('flow')

allMeas = struct2cell(results);
allMeas = cellfun(@(x) mean(x.(usedMeasure),1,'omitnan'),allMeas(1:end-1),'un',0);
allMeas = cat(1,allMeas{:});

ylims = quantile(allMeas,[0 1],'all');
for expType = expTypes
    currentMeasure = results.(expType{1}).(usedMeasure)(:,:,chIdx);
    avgMeas = squeeze(mean(currentMeasure,1,'omitnan'));
    switch usedMeasure
        case {'pec','ppc','coh','icoh','pli','wpli'}
            f = results.(expType{1}).specFreqs;
            pltStr = '-';
        case 'bp'
            f = mean(results.params.freqBands,2);
            pltStr = 'x';
    end
    nexttile;
    hold on
    plot(f,avgMeas,pltStr)
    if ~use_matched_pairs
        plot(f,mean(avgMeas,2),[pltStr 'k'],'LineWidth',3)
    end
    title([strrep(expType{1},'_',' ') ' - ' usedMeasure]);
    ylim(ylims);
    xlim(results.params.freqLims)
    xlabel('Frequency (Hz)')
    ylabel(usedMeasure);
    axis square
    set(gca,'FontSize',15)
    set(gca,'XScale','log')
    set(gca,'XTick',unique([1 2 5*(round(logspace(0,1.75,10)/5))]))
end

