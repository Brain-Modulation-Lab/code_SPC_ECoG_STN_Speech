function fh = create_SPCSTN_CSFS(cfg, DB);
EvtTypes = cfg.cfg.plot.EvtTypes;
EvtTimes = cfg.EvtTimes;
targets = cfg.targets;
T = cfg.T;
nperms = cfg.nperms;
atlas = cfg.atlas;
min_subj = cfg.min_subj;
min_pairs = cfg.min_pairs;
flag_sort = cfg.flag_sort;
bands = cfg.bands;
% load clusters
SPCDensity = readtable(cfg.map);
% load stn
load(cfg.ref,'atlases'); % manually load definition of DISTAL atlas.

keep_labels = {'STN.','STN_motor','STN_associative','STN_limbic'};
idx_labels = cellfun(@(x) find(contains(atlases.names,x)), keep_labels);

STN=arrayfun(@(x) reducepatch(atlases.roi{x,2}.fv,.5),idx_labels,'UniformOutput',false); % extract the left side


ECoG_atlas = DB.Pairs.Location.(atlas);
ECoG_subject = cellstr(DB.Pairs.Location.subj_id);
[ROI_AtlasLabels,ROI_nPairs]  = getROIs(ECoG_atlas, ECoG_subject, min_subj,min_pairs,flag_sort);
% remove some areas like MTG, pars. O trian and STG plan
ROI_AtlasLabels(end-2:end) = [];
ROI_nPairs(end-2:end) = [];


% find coordinates (Y for CS and Z for FS)

cfg.VarFields = {'ClusterCentroidBand','ClusterS_MNI'};
[ClusterCentroidBand, ClusterS_MNI] = get_SPCClustersProperties(cfg,DB.Clusters);

% average and peaks of SPC for each frequency band
ClusterS_MNI_Band_avg = arrayfun(@(x) mean(ClusterS_MNI(ClusterCentroidBand == x,:),'omitnan'),2:6,'uni',false);
ClusterS_MNI_Band_peak = nan(5,3);

for fi = 1:5
    %tmp = SPCDensity(SPCDensity.Color == fi,:);
    %[~,idx_] = max(tmp.Size,[],'omitnan');
    %ClusterS_MNI_Band_peak(fi,:) = tmp{idx_,{'X','Y','Z'}};
    tmp = SPCDensity{:,fi+4};
    [~,idx_] = max(tmp,[],'omitnan');
    ClusterS_MNI_Band_peak(fi,:) = SPCDensity{idx_,{'X_MNI','Y_MNI','Z_MNI'}};
end


for fi = 2:6
    fprintf("Peak coordinate of %s t-SPC: (X,Y,Z) = [%1.2f, %1.2f, %1.2f] mm \n", bands.name{fi},ClusterS_MNI_Band_peak(fi-1,1),ClusterS_MNI_Band_peak(fi-1,2),ClusterS_MNI_Band_peak(fi-1,3))
    fprintf("Average coordinate of %s t-SPC: (X,Y,Z) = [%1.2f, %1.2f, %1.2f] mm \n", bands.name{fi},ClusterS_MNI_Band_avg{fi-1}(1),ClusterS_MNI_Band_avg{fi-1}(2),ClusterS_MNI_Band_avg{fi-1}(3))
end


%% coordinates plot

CoordBands = cell(1,6);
for fi = 2:6
    CoordBands{fi} = ClusterS_MNI(ClusterCentroidBand == fi,:) - mean(STN{1}.vertices);
end

pairs_test = nchoosek(2:6,2);
p_test = cell(1,3);
diff_test = cell(1,3);
for xi = 1 : 3
    p_test{xi} = nan(size(pairs_test,1),1);
    diff_test{xi} =  nan(size(pairs_test,1),1);
    for pi = 1 :  size(pairs_test,1)
        [p_test{xi}(pi), diff_test{xi}(pi)]  = permutation_diffTest_indsamples(CoordBands{pairs_test(pi,1)}(:,xi), CoordBands{pairs_test(pi,2)}(:,xi), 500);
    end
end

fh{1} = figure('position',[200 200 900 600]);
tiledlayout(5,3)
for fi = 2: 6
    nexttile
    raincloud_plot(CoordBands{fi}(:,1), 'box_on', 0, 'color', [0 0 0], 'alpha', 0.5,...
        'box_dodge', 0, 'box_dodge_amount', 0, 'dot_dodge_amount', 0,...
        'box_col_match', 0);
    box off
    hold on
    xline( median(CoordBands{fi}(:,1)),'r--',sprintf('   %1.2f mm ', median(CoordBands{fi}(:,1))),'LabelOrientation','horizontal','linewidth',1.5)
    xlim([-10 10])

    nexttile
    raincloud_plot(CoordBands{fi}(:,2), 'box_on', 0, 'color', [0 0 0], 'alpha', 0.5,...
        'box_dodge', 0, 'box_dodge_amount', 0, 'dot_dodge_amount', 0,...
        'box_col_match', 0);
    box off
    hold on
    xline( median(CoordBands{fi}(:,2)),'r--',sprintf('   %1.2f mm ', median(CoordBands{fi}(:,2))),'LabelOrientation','horizontal','linewidth',1.5)
    xlim([-10 10])

    nexttile
    raincloud_plot(CoordBands{fi}(:,3), 'box_on', 0, 'color', [0 0 0], 'alpha', 0.5,...
        'box_dodge', 0, 'box_dodge_amount', 0, 'dot_dodge_amount', 0,...
        'box_col_match', 0);
    box off
    hold on
    xline( median(CoordBands{fi}(:,3)),'r--',sprintf('   %1.2f mm ', median(CoordBands{fi}(:,3))),'LabelOrientation','horizontal','linewidth',1.5)
    xlim([-10 10])

end


%%

STN_vx = STN{1}.vertices;
[coeff,score,latent] = pca(zscore(STN_vx));

% how to apply transform (((STN_PCtbl{1,1:3} - mean(STN_vx))./std(STN_vx))*coeff).*std(STN_vx)

%scale_fac_dots = 5;
scale_fac_pc = 5;

STN_PCtbl = table(STN_vx(:,1),STN_vx(:,2),STN_vx(:,3),score(:,1)*std(STN_vx(:,1)),score(:,2)*std(STN_vx(:,2)),score(:,3)*std(STN_vx(:,3)),...
    'VariableNames',{'X','Y','Z','PC1','PC2','PC3'});

fh{2} = figure('renderer','painters','position',[300 300 1600 500]);
tiledlayout(1,4)
nexttile
patch('faces',STN{2}.faces, 'vertices',STN{2}.vertices,"FaceColor",[0.9290 0.6940 0.1250],"edgecolor","none","FaceAlpha",0.4)
patch('faces',STN{3}.faces, 'vertices',STN{3}.vertices,"FaceColor",[0 0.4470 0.7410],"edgecolor","none","FaceAlpha",0.4)
patch('faces',STN{4}.faces, 'vertices',STN{4}.vertices,"FaceColor",[1 1 0],"edgecolor","none","FaceAlpha",0.4)
hold on
for pc_i = 1:3
    plot3(mean(STN_vx(:,1)) + scale_fac_pc*sqrt(latent(pc_i))*repmat(coeff(1,pc_i),1,2).*[-1 1], mean(STN_vx(:,2)) + scale_fac_pc*sqrt(latent(pc_i))*repmat(coeff(2,pc_i),1,2).*[-1 1], mean(STN_vx(:,3)) + scale_fac_pc*sqrt(latent(pc_i))*repmat(coeff(3,pc_i),1,2).*[-1 1],'linewidth',4,'color','r');
    %plot3(mean(ClusterE_MNI_valid(:,1)) + coeff(1,pc_i) + param_mapping, mean(ClusterE_MNI_valid(:,2))+ coeff(2,pc_i) + param_mapping, mean(ClusterE_MNI_valid(:,3)) + coeff(3,pc_i) + param_mapping,'linewidth',2,'color',color_pc(pc_i,:))
end
set(gca,"view",1.0e+02 *[ -1.394536101759164   0.089355260649410])
xlabel('X [mm]')
ylabel('Y [mm]')
zlabel('Z [mm]')
box off

for pc_i = 1:3
    nexttile
    s = scatter3(STN_PCtbl.X,STN_PCtbl.Y,STN_PCtbl.Z,15,STN_PCtbl.(['PC',num2str(pc_i)]),'filled');
    %s.SizeData = PC_tbl.nPoints;
    colormap(linspecer)
    view(-60,45)
    xlabel('X [mm]')
    ylabel('Y [mm]')
    zlabel('Z [mm]')
    cb = colorbar;
    ylabel(cb,sprintf('PC%d score [a.u.]',pc_i))
    box off
    set(gca,"view",1.0e+02 *[ -1.394536101759164   0.089355260649410])
    grid off
end

% at the cluster level
ClusterLoc_PC = table();
ClusterLoc_PC.X_MNI = ClusterS_MNI(:,1);
ClusterLoc_PC.Y_MNI = ClusterS_MNI(:,2);
ClusterLoc_PC.Z_MNI = ClusterS_MNI(:,3);

tempPC = nan(height(ClusterLoc_PC),3);
for ii = 1 : height(ClusterLoc_PC)
    tt = ClusterLoc_PC{ii,{'X_MNI','Y_MNI','Z_MNI'}};
    tempPC(ii,:) = (tt - mean(STN_vx))./std(STN_vx)*coeff.*std(STN_vx);

end
ClusterLoc_PC.PC1 = tempPC(:,1);
ClusterLoc_PC.PC2 = tempPC(:,2);
ClusterLoc_PC.PC3 = tempPC(:,3);

% at the density level

tempPC = nan(height(SPCDensity),3);

for ii = 1 : height(SPCDensity)
    tt = SPCDensity{ii,{'X_MNI','Y_MNI','Z_MNI'}};
    tempPC(ii,:) = (tt - mean(STN_vx))./std(STN_vx)*coeff.*std(STN_vx);

end
SPCDensity.PC1 = tempPC(:,1);
SPCDensity.PC2 = tempPC(:,2);
SPCDensity.PC3 = tempPC(:,3);

% all pairs
PairLoc_PC = table();
PairLoc_PC.X_MNI = DB.Pairs.Location{:,{'S_MNI_X'}};
PairLoc_PC.Y_MNI = DB.Pairs.Location{:,{'S_MNI_Y'}};
PairLoc_PC.Z_MNI = DB.Pairs.Location{:,{'S_MNI_Z'}};


tempPC = nan(height(PairLoc_PC),3);
for ii = 1 : height(PairLoc_PC)
    tt = PairLoc_PC{ii,{'X_MNI','Y_MNI','Z_MNI'}};
    tempPC(ii,:) = (tt - mean(STN_vx))./std(STN_vx)*coeff.*std(STN_vx);
end
PairLoc_PC.PC1 = tempPC(:,1);
PairLoc_PC.PC2 = tempPC(:,2);
PairLoc_PC.PC3 = tempPC(:,3);
% ##########################3
%   this is wrong!!!!!!!
% idx_PC2clust= dsearchn(STN_vx,delaunayn(STN_vx),ClusterLoc_PC{:,{'X_MNI','Y_MNI','Z_MNI'}});
%
% ClusterLoc_PC.PC1 = STN_PCtbl.PC1(idx_PC2clust);
% ClusterLoc_PC.PC2 = STN_PCtbl.PC2(idx_PC2clust);
% ClusterLoc_PC.PC3 = STN_PCtbl.PC3(idx_PC2clust);
% #############################


%% PC coordinates plot

CLUSTER_PC = [ClusterLoc_PC.PC1 ClusterLoc_PC.PC2 ClusterLoc_PC.PC3];
CoordBands = cell(1,6);
for fi = 2:6
    CoordBands{fi} = CLUSTER_PC(ClusterCentroidBand == fi,:);
end


centroid_bands = [-13.54,-15.38, -6.3; -13.23 -15.26 -8.07; -13.32 -14.16 -7.64];
centroid_bands_PC = nan(3,3);
for cc = 1 : 3
    centroid_bands_PC(cc,:) =  (centroid_bands(cc,:) - mean(STN_vx))./std(STN_vx)*coeff.*std(STN_vx);
end

% only spatial information
%[unique_spc1,~,idx_unique1] = unique(SPCDensity{:,{'PC1','PC2'}},'rows');


% find covnex hull (boundary of stn) in PC
STN_convexh = convhull(STN_vx);
STN_convexh_boundaryppc = nan(size(STN_convexh,1),3);
for cc = 1 : size(STN_convexh,1)
    STN_convexh_boundaryppc(cc,:) = (STN_vx(STN_convexh(cc,1),:) - mean(STN_vx))./std(STN_vx)*coeff.*std(STN_vx);
end



fh{3} = figure('Position',[200 200 700 700],'renderer','painter');
%bubblepie(SPCDensity.PC1,SPCDensity.PC2,(SPCDensity.theta + SPCDensity.alpha)/100, SPCDensity{:,{'theta','alpha'}}./sum(SPCDensity{:,{'theta','alpha'}},2,'omitnan'),[],'PC1','PC2', 0,{'#D72027','#FAAE61'});
bubblepie(SPCDensity.PC1,SPCDensity.PC2,SPCDensity.broadband/125,[SPCDensity{:,{'theta','alpha','beta'}} (SPCDensity.lowGamma + SPCDensity.highGamma)]./SPCDensity{:,{'broadband'}},[],'PC1','PC2', 0,true,{'#D72027','#FAAE61','#FFCC5C','#ECECEC'});
hold on
kh = convhull(STN_convexh_boundaryppc(:,1),STN_convexh_boundaryppc(:,2));
plot3(STN_convexh_boundaryppc(kh,1),STN_convexh_boundaryppc(kh,2),ones(numel(kh),1),'k')
dscatter(PairLoc_PC.PC1,PairLoc_PC.PC2,'BINS',[60 60],'plottype','surf');
xlabel('PC1 [mm]')
ylabel('PC2 [mm]')
yline(0,'--k')
xline(0,'--k')
colormap(flipud(colormap('gray')));
%colorbar;
xlim([-5 2.5])
ylim([-6.5 4.5])
clim([0.000000004187411  0.001888804770069]);
%crange = get(gca,'CLim')
% thr_c = crange(2)*100/height(PairLoc_PC);
% clim([thr_c crange(2)]);
view(2)
%saveFigures(gcf,'STN_PC1vsPC2_thetaalpha_beta');






fh{4} = figure('Position',[200 200 700 700],'renderer','painter');
%bubblepie(SPCDensity.PC1,SPCDensity.PC3,(SPCDensity.theta + SPCDensity.alpha)/100, SPCDensity{:,{'theta','alpha'}}./sum(SPCDensity{:,{'theta','alpha'}},2,'omitnan'),[],'PC1','PC2', 0,{'#D72027','#FAAE61'});
bubblepie(SPCDensity.PC1,SPCDensity.PC3,SPCDensity.broadband/125,[SPCDensity{:,{'theta','alpha','beta'}} ( SPCDensity.lowGamma + SPCDensity.highGamma)]./SPCDensity{:,{'broadband'}},[],'PC1','PC3', 0,true,{'#D72027','#FAAE61','#FFCC5C','#ECECEC'});

hold on
kh = convhull(STN_convexh_boundaryppc(:,1),STN_convexh_boundaryppc(:,3));
plot3(STN_convexh_boundaryppc(kh,1),STN_convexh_boundaryppc(kh,3),ones(numel(kh),1),'k')
dscatter(PairLoc_PC.PC1,PairLoc_PC.PC3,'BINS',[60 60],'plottype','surf');
xlabel('PC1 [mm]')
ylabel('PC3 [mm]')
yline(0,'--k')
xline(0,'--k')
colormap(flipud(colormap('gray')));
%colorbar;
% crange = get(gca,'CLim');
% thr_c = crange(2)*100/height(PairLoc_PC);
% clim([thr_c crange(2)]);
xlim([-5 2.5])
ylim([-2.6 2.6])
view(2)
clim([0.000000004187411  0.001888804770069]);
%saveFigures(gcf,'STN_PC1vsPC3_thetaalpha_beta');





fh{5} = figure('Position',[200 200 700 700],'renderer','painter');
%bubblepie(SPCDensity.PC2,SPCDensity.PC3,(SPCDensity.theta + SPCDensity.alpha)/100, SPCDensity{:,{'theta','alpha'}}./sum(SPCDensity{:,{'theta','alpha'}},2,'omitnan'),[],'PC1','PC2', 0,true,{'#D72027','#FAAE61'});
bubblepie(SPCDensity.PC2,SPCDensity.PC3,SPCDensity.broadband/125,[SPCDensity{:,{'theta','alpha','beta'}} ( SPCDensity.lowGamma + SPCDensity.highGamma)]./SPCDensity{:,{'broadband'}},[],'PC2','PC3', 0,true,{'#D72027','#FAAE61','#FFCC5C','#ECECEC'});

hold on
kh = convhull(STN_convexh_boundaryppc(:,2),STN_convexh_boundaryppc(:,3));
plot3(STN_convexh_boundaryppc(kh,2),STN_convexh_boundaryppc(kh,3),ones(numel(kh),1),'k')

dscatter(PairLoc_PC.PC2,PairLoc_PC.PC3,'BINS',[60 60],'plottype','surf');
xlabel('PC2 [mm]')
ylabel('PC3 [mm]')
xlim([-5 5])
ylim([-2.6 2.6])
yline(0,'--k')
xline(0,'--k')
colormap(flipud(colormap('gray')));
%colorbar;
% crange = get(gca,'CLim');
% thr_c = crange(2)*100/height(PairLoc_PC);
% clim([thr_c crange(2)]);
xlim([-6.5 4.5])
ylim([-2.6 2.6])
view(2)
clim([0.000000004187411  0.001888804770069]);





pairs_test = nchoosek(2:6,2);
p_test = cell(1,3);
diff_test = cell(1,3);
for xi = 1 : 3
    p_test{xi} = nan(size(pairs_test,1),1);
    diff_test{xi} =  nan(size(pairs_test,1),1);
    for pi = 1 :  size(pairs_test,1)
        [p_test{xi}(pi), diff_test{xi}(pi)]  = permutation_diffTest_indsamples(CoordBands{pairs_test(pi,1)}(:,xi), CoordBands{pairs_test(pi,2)}(:,xi), 500);
    end
end


fh{6} = figure('position',[200 200 900 600]);
tiledlayout(5,3)
for fi = 2: 6
    nexttile
    raincloud_plot(CoordBands{fi}(:,1), 'box_on', 0, 'color', [0 0 0], 'alpha', 0.5,...
        'box_dodge', 0, 'box_dodge_amount', 0, 'dot_dodge_amount', 0,...
        'box_col_match', 0);
    box off
    hold on
    xline( median(CoordBands{fi}(:,1)),'r--',sprintf('   %1.2f mm ', median(CoordBands{fi}(:,1))),'LabelOrientation','horizontal','linewidth',1.5)
    xlim([-10 10])

    nexttile
    raincloud_plot(CoordBands{fi}(:,2), 'box_on', 0, 'color', [0 0 0], 'alpha', 0.5,...
        'box_dodge', 0, 'box_dodge_amount', 0, 'dot_dodge_amount', 0,...
        'box_col_match', 0);
    box off
    hold on
    xline( median(CoordBands{fi}(:,2)),'r--',sprintf('   %1.2f mm ', median(CoordBands{fi}(:,2))),'LabelOrientation','horizontal','linewidth',1.5)
    xlim([-10 10])

    nexttile
    raincloud_plot(CoordBands{fi}(:,3), 'box_on', 0, 'color', [0 0 0], 'alpha', 0.5,...
        'box_dodge', 0, 'box_dodge_amount', 0, 'dot_dodge_amount', 0,...
        'box_col_match', 0);
    box off
    hold on
    xline( median(CoordBands{fi}(:,3)),'r--',sprintf('   %1.2f mm ', median(CoordBands{fi}(:,3))),'LabelOrientation','horizontal','linewidth',1.5)
    xlim([-10 10])

end

%








% #####################################3
% this is wrong!!!!!!!!
% idx_PC2pair = dsearchn(STN_vx,delaunayn(STN_vx),DB.Pairs.Location{:,{'S_MNI_X','S_MNI_Y','S_MNI_Z'}});
% PairLoc_PC.PC1 = STN_PCtbl.PC1(idx_PC2pair);
% PairLoc_PC.PC2 = STN_PCtbl.PC2(idx_PC2pair);
% PairLoc_PC.PC3 = STN_PCtbl.PC3(idx_PC2pair);
% #####################################3






% check focality
fh{7} = figure('position',[200 200 1600 400]);
tiledlayout(1,5)
for fi = 2:6
    nexttile
    coords = ClusterS_MNI(ClusterCentroidBand == fi,:);
    true_dist_clust = mean(vecnorm(coords - mean(coords,'omitnan'),2,2),'omitnan'); % mm
    perm_dist_clust = nan(1,nperms);
    for perm_i = 1 : nperms
        perm_loc = ClusterS_MNI(randi(height(ClusterS_MNI),1,size(coords,1)),:);
        perm_dist_clust(perm_i) = mean(vecnorm(perm_loc - mean(perm_loc,'omitnan'),2,2),'omitnan'); % mm
    end
    histogram(perm_dist_clust,round(sqrt(nperms)),'FaceColor',[.8 ,.8 ,.8])
    hold on
    box off
    xlabel('Distance to centre [mm]')
    ylabel(' # observations ')
    xlim([0.8 3])
    if true_dist_clust > min(perm_dist_clust)
        p_val = (sum((perm_dist_clust - mean(perm_dist_clust)) <= (true_dist_clust - mean(perm_dist_clust)))+1) / (nperms+1); %V2 from Ernst, Permutation Methods: A Basis for Exact Inference, 2004
    else
        p_val = 1/nperms;
    end
    xline(true_dist_clust,'label',{sprintf('p = %1.3f',p_val)},'Color','r','linewidth',1.4,'LabelOrientation','horizontal','FontSize',15)
end
% 3d plot
%%


target = {'STN','GPi'};
faceColors = {'#F7B05B','#B2E4DB'};
faceAlphas = [.6 .6];

atlas_names = cellfun(@(x) x.name, atlases.roi, 'UniformOutput', false);
[~,atlas_idx] = ismember(target, atlas_names);
assert(any(~isempty(atlas_idx)),"No target in this atlas. Change the Atlas!");

faces = arrayfun(@(x) atlases.roi{x,2}.fv.faces,atlas_idx,'UniformOutput',false);
vertices = arrayfun(@(x) atlases.roi{x,2}.fv.vertices,atlas_idx,'UniformOutput',false);

% here I need to quantify how many neurons and create an electroed table


for fi = 1:3
    fh{7 + fi} = figure('position',[200 200 1600 400]);
    tiledlayout(1,2)
    electrode = table();
    electrode.x = [SPCDensity.X_MNI];
    electrode.y = SPCDensity.Y_MNI;
    electrode.z = SPCDensity.Z_MNI;
    electrode.radius = log(SPCDensity{:,fi + 4} + 1)*1.01 + 0.1; % scale down for visualization
    electrode.color = repmat(hex2rgb(bands.color(fi + 1)),height(SPCDensity),1);

    nexttile
    cfg = [];
    cfg.h_ax = gca;
    cfg.view = [-20,15];
    %title({[SUBJECT], ' native-space ECoG reconstruction'})
    hold on
    for target_i = 1 : numel(target)
        cfg.surface_facecolor = hex2rgb(faceColors{target_i});
        cfg.surface_facealpha = faceAlphas(target_i);

        bml_plot3d_surface(cfg, vertices{target_i},faces{target_i});
    end

    hold on

    cfg.h_ax = gca;
    cfg.annotate = false;
    bml_plot3d_points(cfg, electrode);
    nexttile
    cfg = [];
    cfg.h_ax = gca;
    cfg.view = [20,15];
    %title({[SUBJECT], ' native-space ECoG reconstruction'})
    hold on
    for target_i = 1 : numel(target)
        cfg.surface_facecolor = hex2rgb(faceColors{target_i});
        cfg.surface_facealpha = faceAlphas(target_i);

        bml_plot3d_surface(cfg, vertices{target_i},faces{target_i});
    end

    hold on

    cfg.h_ax = gca;
    cfg.annotate = false;
    bml_plot3d_points(cfg, electrode);
end







