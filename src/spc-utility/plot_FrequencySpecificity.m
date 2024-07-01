function fh = plot_FrequencySpecificity(cfg,DB)

% get cfg fields
bands = cfg.bands;
freq_vec = mean(cfg.freqsOfInt,2);
nPerms = cfg.nPerms;

% get PPCz power at the single-pair level
nPairs = numel(DB.Pairs.PPCzmap);
thetaPPCz = nan(1, nPairs);
alphaPPCz = nan(1, nPairs);
betaPPCz = nan(1, nPairs);
for pair_i = 1 : nPairs
    thetaPPCz(pair_i) = mean(mean(DB.Pairs.PPCzmap{pair_i}(freq_vec <= bands.fends(2) & freq_vec >= bands.fstarts(2),:),1,'omitnan'),2,'omitnan');
    alphaPPCz(pair_i) = mean(mean(DB.Pairs.PPCzmap{pair_i}(freq_vec <= bands.fends(3) & freq_vec >= bands.fstarts(3),:),1,'omitnan'),2,'omitnan');
    betaPPCz(pair_i) = mean(mean(DB.Pairs.PPCzmap{pair_i}(freq_vec <= bands.fends(4) & freq_vec >= bands.fstarts(4),:),1,'omitnan'),2,'omitnan');
end

% get t-SPC frequency dendity at unit level 
[SPCbandsunits_Density,BandSpecifity_btsp] = get_SingleUnitDensity(DB, nPerms);
nUnits = height(SPCbandsunits_Density);

fh{1} = figure;
tiledlayout(2,3)
nexttile(1,[1 3])
jitter = 1 + 0.2*randn(1,nUnits);
symbol = {'v','x','^','o'};
type = [-1 0 1 2];
scatter(SPCbandsunits_Density.band_spec,jitter,35,'k','filled')
ylim([-1.75 2.75])
xlim([-.2 1.2])
xline(prctile(BandSpecifity_btsp,5))
xline(prctile(BandSpecifity_btsp,50))
xline(prctile(BandSpecifity_btsp,95))

mean(SPCbandsunits_Density.band_spec,'omitnan')
std(SPCbandsunits_Density.band_spec,'omitnan')/sqrt(nUnits)
xlabel('Frequency-Specificity')
yticks([])

nexttile(4)
idx_ = dsearchn(SPCbandsunits_Density.band_spec, 0.0426713);
color = [1 0  0 ; 0 1 0; 0 0 1];
p1 = pie(SPCbandsunits_Density{idx_,{'theta','alpha','beta'}});
p1(1).FaceColor = color(1,:);
p1(3).FaceColor = color(2,:);
p1(5).FaceColor = color(3,:);

nexttile(5)
idx_ = 87;%dsearchn(SPCbandsunits_Density.band_spec, 1);
color = [1 0  0 ; 0 1 0; 0 0 1; 1 1 0; 0 1 1];
p1 = pie(SPCbandsunits_Density{idx_,{'theta','alpha','beta'}});
p1(1).FaceColor = color(1,:);
p1(3).FaceColor = color(2,:);
p1(5).FaceColor = color(3,:);

nexttile(6)
idx_ = dsearchn(SPCbandsunits_Density.band_spec, 0.65);
color = [1 0  0 ; 0 1 0; 0 0 1; 1 1 0; 0 1 1];
p1 = pie(SPCbandsunits_Density{idx_,{'theta','alpha','beta'}});
p1(1).FaceColor = color(1,:);
p1(3).FaceColor = color(2,:);
p1(5).FaceColor = color(3,:);



fh{2} = figure('Position',[200 200 1200 350]);
tl = tiledlayout(1,3);
nexttile(tl,1);
dscatter(thetaPPCz',alphaPPCz')
hold on
dscatter(thetaPPCz',alphaPPCz','plottype','contour')
colormap(linspecer)
%colorbar
xlim([-1 8])
ylim([-.3 1.7])
xlabel('\theta PPC [z-score]')
ylabel('\alpha PPC [z-score]')


nexttile(tl,2)
dscatter(thetaPPCz',betaPPCz')
hold on
dscatter(thetaPPCz',betaPPCz','plottype','contour')
colormap(linspecer)
%colorbar
xlim([-1 8])
ylim([-.3 1.7])
xlabel('\theta PPC [z-score]')
ylabel('\beta PPC [z-score]')


nexttile(tl,3)
dscatter(alphaPPCz',betaPPCz')
hold on
dscatter(alphaPPCz',betaPPCz','plottype','contour')
colormap(linspecer)
%colorbar
xlim([-1 8])
ylim([-.3 1.7])
xlabel('\alpha PPC [z-score]')
ylabel('\beta PPC [z-score]')