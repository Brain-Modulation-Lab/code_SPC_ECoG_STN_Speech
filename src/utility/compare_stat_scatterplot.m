function varargout = compare_stat_scatterplot(x,y,type)

%remove nans
isnan_x = isnan(x);
isnan_y = isnan(y);


switch type
    case 'corr'
        x = x(~isnan_x & ~isnan_y);
        y = y(~isnan_x & ~isnan_y);
        [Rho,pRho] = permucorr(x(:),y(:),"verbose",0);
        varargout{1} = Rho;
        varargout{2} = pRho;
    case 'corr_sp'
        x = x(~isnan_x & ~isnan_y);
        y = y(~isnan_x & ~isnan_y);
        [Rho,pRho] = permucorr(x(:),y(:),"verbose",0,'type','spearman');
        varargout{1} = Rho;
        varargout{2} = pRho;        

    case 'ttest_dep'
        x = x(~isnan_x & ~isnan_y);
        y = y(~isnan_x & ~isnan_y);
        [T,pT] = permuttest(x(:),y(:),"verbose",0);
        varargout{1} = T;
        varargout{2} = pT;

    case 'ttest_ind'
        x = x(~isnan_x);
        y = y(~isnan_y);
        [T,pT] = permuttest2(x(:),y(:),"verbose",0);
        varargout{1} = T;
        varargout{2} = pT;

end

