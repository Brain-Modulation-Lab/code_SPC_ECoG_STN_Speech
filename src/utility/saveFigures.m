function saveFigures(fh, figname, ext)

if ~exist('ext', 'var')
    saveas(fh,[figname, '.fig'])
    saveas(fh,[figname, '.png'])
    saveas(fh,[figname, '.pdf'])
    saveas(fh,[figname, '.svg'])
    exportgraphics(fh,[figname, '.eps'],'ContentType','vector')

    
else
    if ischar(ext)
        tmp = ext;
        ext = cell(1);
        ext{1} = tmp;
    end
    for ext_i = 1 : numel(ext)
        if strcmpi(ext{ext_i},'.eps')
            exportgraphics(fh,[figname, '.eps'],'ContentType','vector')
        else
            saveas(fh,[figname, ext{ext_i}])
        end
    end
end
