function ax = scatter_jitter_paired(data,labels,dispersion,spacing,color)


if ismatrix(color)
    color = mat2cell(color,[1 1])';
end
assert(numel(data{1}) == numel(data{2}),'Observations must be paired!')

offset = 1;

hold on
%data = reshape([ROI.fr{:}],[],1);

for win_i = 1 : numel(labels)
    jitter = dispersion*randn(1,numel(data{win_i}));
    scatter(jitter + spacing(win_i),data{win_i},25,'MarkerEdgeColor',color{win_i},'LineWidth',2, ...
        'MarkerFaceColor',color{win_i},'MarkerFaceAlpha',.7)
    line(spacing(win_i) + [-.2 .2],median(data{win_i})*ones(1,2),'color','k','linewidth',3)
end

% creating connecting lines
line_ext1 = spacing(1) + diff(spacing)/5;
line_ext2 = spacing(2) - diff(spacing)/5;

line(repmat([line_ext1; line_ext2],1,numel(data{1})),[data{1} data{2}]','color',[.4 .4 .4])
xticks(spacing)
xlim([spacing(1) - (2*dispersion + offset) spacing(end) + (2*dispersion + offset)])
%ylabel('Firing rate [z-score]')
xticklabels(labels)

ax = gca;

end
