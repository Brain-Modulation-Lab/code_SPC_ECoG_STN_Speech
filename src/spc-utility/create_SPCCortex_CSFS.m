function fh = create_SPCCortex_CSFS(cfg, DB);

bands = cfg.bands;
nperms = cfg.nperms;
atlas = cfg.atlas;
min_subj = cfg.min_subj;
min_pairs = cfg.min_pairs;
flag_sort = cfg.flag_sort;
% load clusters

SPCDensity = readtable(cfg.map);
% load cortex
load(cfg.ref);


% find coordinates (Y for CS and Z for FS)

cfg.VarFields = {'ClusterCentroid','ClusterCentroidBand','ClusterE_MNI','ClusterEAtlasLabel','ClusterDur'};
[ClusterCentroid,ClusterCentroidBand, ClusterE_MNI,ClusterEAtlasLabel,ClusterDur] = get_SPCClustersProperties(cfg,DB.Clusters);

% average and peaks of SPC for each frequency band
ClusterS_MNI_Band_avg = arrayfun(@(x) mean(ClusterE_MNI(ClusterCentroidBand == x,:),'omitnan'),2:6,'uni',false);
ClusterS_MNI_Band_peak = nan(5,3);

for fi = 1:5
    tmp = SPCDensity(SPCDensity.Color == fi,:);
    [~,idx_] = max(tmp.Size,[],'omitnan');
    ClusterS_MNI_Band_peak(fi,:) = tmp{idx_,{'X','Y','Z'}};
end

fh{1} = figure; % 2d viz
for fi = 1:5
    idx_ = SPCDensity.Color == fi;
    hold on
    scatter(SPCDensity.Y(idx_), SPCDensity.Z(idx_),SPCDensity.Size(idx_)*500, ...
        'MarkerFaceColor',hex2rgb(cfg.bands.color(fi+1)), ...
        'MarkerEdgeColor',hex2rgb(cfg.bands.color(fi+1)), ...
        'MarkerFaceAlpha',.3)
    scatter(ClusterS_MNI_Band_avg{fi}(2),ClusterS_MNI_Band_avg{fi}(3), 400, ...
        'MarkerFaceColor',hex2rgb(cfg.bands.color(fi+1)), ...
        'MarkerEdgeColor',hex2rgb(cfg.bands.color(fi+1)), ...
        'MarkerFaceAlpha',1,'Marker','^','linewidth',3)
end

 set ( gca, 'xdir', 'reverse' )
 xlabel(' Y [mm]')
ylabel(' Z [mm]')

%atlas_region = compute_MNIcoords(ClusterS_MNI_Band_peak,BS1);



% check density, freq, duration and duration
% get ROIs
ECoG_atlas = DB.Pairs.Location.(atlas);
ECoG_subject = cellstr(DB.Pairs.Location.subj_id);
[ROI_AtlasLabels,ROI_nPairs]  = getROIs(ECoG_atlas, ECoG_subject, min_subj,min_pairs,flag_sort);
% remove some areas like MTG, pars. O trian and STG plan
ROI_AtlasLabels(end-2:end) = [];
ROI_nPairs(end-2:end) = [];
% density
ROI_nClusters = cellfun(@(x)sum(strcmpi(ClusterEAtlasLabel,x)),ROI_AtlasLabels);
ROI_SPCdensity = ROI_nClusters./ROI_nPairs*100;


% need build permutation
ROI_SPCdensity_perm = cell(1,numel(ROI_AtlasLabels));
ROI_SPCdensity_perm_5pth = nan(1,numel(ROI_AtlasLabels));
ROI_SPCdensity_perm_50pth = nan(1,numel(ROI_AtlasLabels));
ROI_SPCdensity_perm_95pth = nan(1,numel(ROI_AtlasLabels));

for roi_i = 1 : numel(ROI_AtlasLabels)
    for perm_i = 1 : nperms
        tmp = DB.Pairs.Location.nSign(randperm(height(DB.Pairs.Location), ROI_nPairs(roi_i)));
        ROI_SPCdensity_perm{roi_i}(perm_i) = mean(tmp)*100;
    end
    ROI_SPCdensity_perm_5pth(roi_i) = prctile(ROI_SPCdensity_perm{roi_i},5);
    ROI_SPCdensity_perm_50pth(roi_i) = prctile(ROI_SPCdensity_perm{roi_i},50);
    ROI_SPCdensity_perm_95pth(roi_i) = prctile(ROI_SPCdensity_perm{roi_i},95);
    
end

fh{2} = figure('position',[200 200 400 400]); % density
bar(1:numel(ROI_AtlasLabels),ROI_SPCdensity)
xticks(1:numel(ROI_AtlasLabels))
hold on
plot(1:numel(ROI_AtlasLabels),ROI_SPCdensity_perm_5pth,'o')
plot(1:numel(ROI_AtlasLabels),ROI_SPCdensity_perm_50pth,'o')
plot(1:numel(ROI_AtlasLabels),ROI_SPCdensity_perm_95pth,'o')

xticklabels(ROI_AtlasLabels)
ylabel('t-SPC density')
box off

% duration
fh{3} = figure('position',[200 200 400 1000]);
ROI_ClustDur = cellfun(@(x) ClusterDur(strcmpi(ClusterEAtlasLabel,x)),ROI_AtlasLabels,'uni',false);
tiledlayout(numel(ROI_AtlasLabels),1)
for roi_i = 1 : numel(ROI_AtlasLabels)
    nexttile
    raincloud_plot(ROI_ClustDur{roi_i}, 'box_on', 0, 'color', [0 0 0], 'alpha', 0.5,...
        'box_dodge', 0, 'box_dodge_amount', 0, 'dot_dodge_amount', 0,...
        'box_col_match', 0);
    box off
    hold on
    xline( median(ROI_ClustDur{roi_i}),'r--',sprintf('   %1.2f s ', median(ROI_ClustDur{roi_i})),'LabelOrientation','horizontal','linewidth',1.5)

xlim([0 1])
ylim([0 4.8])
title(ROI_AtlasLabels{roi_i})
end

% freq
fh{4} = figure('position',[200 200 400 1000]);
ROI_ClustFreq = cellfun(@(x) ClusterCentroid(strcmpi(ClusterEAtlasLabel,x),2),ROI_AtlasLabels,'uni',false);
tiledlayout(numel(ROI_AtlasLabels),1)
for roi_i = 1 : numel(ROI_AtlasLabels)
    nexttile
    raincloud_plot(ROI_ClustFreq{roi_i}, 'box_on', 0, 'color', [0 0 0], 'alpha', 0.5,...
        'box_dodge', 0, 'box_dodge_amount', 0, 'dot_dodge_amount', 0,...
        'box_col_match', 0);
    box off
    hold on
    xline( median(ROI_ClustFreq{roi_i}),'r--',sprintf('   %1.2f Hz ', median(ROI_ClustFreq{roi_i})),'LabelOrientation','horizontal','linewidth',1.5)

    set(gca,'xscale','log')
     ylim([0 0.055])
    xlim([4 150])
    title(ROI_AtlasLabels{roi_i})
end



% check focality

fh{5} = figure('position',[200 200 1800 300]);
tiledlayout(1,5)

for fi = 1:5
    coords = ClusterE_MNI(ClusterCentroidBand == fi+1,:);
    true_dist_clust = mean(vecnorm(coords - mean(coords,'omitnan'),2,2),'omitnan'); % mm
    perm_dist_clust = nan(1,nperms);
    for perm_i = 1 : nperms
        perm_loc = ClusterE_MNI(randi(height(ClusterE_MNI),1,size(coords,1)),:);
        perm_dist_clust(perm_i) = mean(vecnorm(perm_loc - mean(perm_loc,'omitnan'),2,2),'omitnan'); % mm
    end
    nexttile
    histogram(perm_dist_clust,round(sqrt(nperms)),'FaceColor',[.8 ,.8 ,.8])
    hold on
    box off
    xlabel('Distance to centre [mm]')
    ylabel(' # observations ')
    xlim([0 30])
    if true_dist_clust > min(perm_dist_clust)
        p_val = (sum((perm_dist_clust - mean(perm_dist_clust)) <= (true_dist_clust - mean(perm_dist_clust)))+1) / (nperms+1); %V2 from Ernst, Permutation Methods: A Basis for Exact Inference, 2004
    else
        p_val = 1/nperms;
    end
    xline(true_dist_clust,'label',{sprintf('p = %1.3f',p_val)},'Color','r','linewidth',1.4,'LabelOrientation','horizontal','FontSize',15)
end

% check unit_Type contribution to roi

ROI_UnitFRmod_density = nan(4,numel(ROI_AtlasLabels));
ROI_UnitFRmod_Pairs = arrayfun(@(x)DB.Pairs.Location(DB.Pairs.Location.S_typeFRmod == x,:),[-1 0 1 2],'uni',false);
for unit_i = 1 : 4
    tmp = ROI_UnitFRmod_Pairs{unit_i};
    for roi_i = 1 : numel(ROI_AtlasLabels)
    ROI_UnitFRmod_density(unit_i,roi_i) = 100*mean(tmp.nSign(strcmpi(tmp.E_atlas_label_Destrieux,ROI_AtlasLabels{roi_i})));
    end
end


ROI_UnitFRmod_SPCdensity_perm = cell(4,numel(ROI_AtlasLabels));
ROI_UnitFRmod_SPCdensity_perm_5pth = nan(4,numel(ROI_AtlasLabels));
ROI_UnitFRmod_SPCdensity_perm_50pth = nan(4,numel(ROI_AtlasLabels));
ROI_UnitFRmod_SPCdensity_perm_95pth = nan(4,numel(ROI_AtlasLabels));

for unit_i = 1:4
    for roi_i = 1 : numel(ROI_AtlasLabels)
        tmp = ROI_UnitFRmod_Pairs{unit_i};

        for perm_i = 1 : nperms
            tmp2 = tmp.nSign(randperm(height(tmp), sum(strcmpi(tmp.E_atlas_label_Destrieux,ROI_AtlasLabels{roi_i}))));
            ROI_UnitFRmod_SPCdensity_perm{unit_i,roi_i}(perm_i) = mean(tmp2)*100;
        end
        ROI_UnitFRmod_SPCdensity_perm_5pth(unit_i,roi_i) = prctile(ROI_UnitFRmod_SPCdensity_perm{unit_i,roi_i},5);
        ROI_UnitFRmod_SPCdensity_perm_50pth(unit_i,roi_i) = prctile(ROI_UnitFRmod_SPCdensity_perm{unit_i,roi_i},50);
        ROI_UnitFRmod_SPCdensity_perm_95pth(unit_i,roi_i) = prctile(ROI_UnitFRmod_SPCdensity_perm{unit_i,roi_i},95);

    end
end
% 
% fh{6} = figure
% plot(ROI_UnitFRmod_density')
% xticks(1:numel(ROI_AtlasLabels))
% xticklabels(ROI_AtlasLabels)
% box off
% ylabel('t-SPC density')
% legend({'D','N','I','M'})



fh{6} = figure('position',[200 200 1200 400]);
tiledlayout(1,4)
for unit_i = 1:4
    nexttile
    bar(1:numel(ROI_AtlasLabels),ROI_UnitFRmod_density(unit_i,:))
    xticks(1:numel(ROI_AtlasLabels))
    hold on
    plot(1:numel(ROI_AtlasLabels),ROI_UnitFRmod_SPCdensity_perm_5pth(unit_i,:),'o')
    plot(1:numel(ROI_AtlasLabels),ROI_UnitFRmod_SPCdensity_perm_50pth(unit_i,:),'o')
    plot(1:numel(ROI_AtlasLabels),ROI_UnitFRmod_SPCdensity_perm_95pth(unit_i,:),'o')

    xticklabels(ROI_AtlasLabels)
    ylabel('t-SPC density')
    box off
    ylim([0 45])
    
end



%% check divided by frequency bands

ROI_FrequencyBand_density = nan(6,numel(ROI_AtlasLabels));
%ROI_UnitFRmod_Pairs = arrayfun(@(x)DB.Pairs.Location(DB.Pairs.Location.S_typeFRmod == x,:),[-1 0 1 2],'uni',false);

for roi_i = 1 : numel(ROI_AtlasLabels)
    for fi = 2:6
        ROI_FrequencyBand_density(fi,roi_i) = 100*sum(ClusterCentroidBand == fi & strcmpi(ClusterEAtlasLabel',ROI_AtlasLabels{roi_i}))/ROI_nPairs(roi_i);
    end
end



%need build permutation
ROI_FrequencyBand_SPCdensity_perm = cell(6,numel(ROI_AtlasLabels));
ROI_FrequencyBand_SPCdensity_perm_5pth = nan(6,numel(ROI_AtlasLabels));
ROI_FrequencyBand_SPCdensity_perm_50pth = nan(6,numel(ROI_AtlasLabels));
ROI_FrequencyBand_SPCdensity_perm_95pth = nan(6,numel(ROI_AtlasLabels));



for roi_i = 1 : numel(ROI_AtlasLabels)
    tmp = ClusterCentroidBand(strcmpi(ClusterEAtlasLabel,ROI_AtlasLabels{roi_i}));
    for fi = 2:6
        for perm_i = 1 : nperms
            tmp2 = ClusterCentroidBand(randperm(numel(ClusterCentroidBand), numel(tmp)));

            ROI_FrequencyBand_SPCdensity_perm{fi,roi_i}(perm_i) = 100*sum(tmp2 == fi)/ROI_nPairs(roi_i);
        end
        ROI_FrequencyBand_SPCdensity_perm_5pth(fi,roi_i) = prctile(ROI_FrequencyBand_SPCdensity_perm{fi,roi_i},5);
        ROI_FrequencyBand_SPCdensity_perm_50pth(fi,roi_i) = prctile(ROI_FrequencyBand_SPCdensity_perm{fi,roi_i},50);
        ROI_FrequencyBand_SPCdensity_perm_95pth(fi,roi_i) = prctile(ROI_FrequencyBand_SPCdensity_perm{fi,roi_i},95);
    end
end






%
fh{7} = figure('position',[200 200 1600 400]);
tiledlayout(1,5)
for fi = 2:6
    nexttile
    bar(1:numel(ROI_AtlasLabels),ROI_FrequencyBand_density(fi,:))
    xticks(1:numel(ROI_AtlasLabels))
    hold on
     plot(1:numel(ROI_AtlasLabels),ROI_FrequencyBand_SPCdensity_perm_5pth(fi,:),'o')
     plot(1:numel(ROI_AtlasLabels),ROI_FrequencyBand_SPCdensity_perm_50pth(fi,:),'o')
   plot(1:numel(ROI_AtlasLabels),ROI_FrequencyBand_SPCdensity_perm_95pth(fi,:),'o')

    xticklabels(ROI_AtlasLabels)
    ylabel('t-SPC density')
    box off
    ylim([0 15.2])
    
end

% 3d figure 


% cortex 3d
cortex = reducepatch(BS1.Faces, BS1.Vertices, 0.3); 


faces = cortex.faces;
vertices = cortex.vertices;

electrode = table();
electrode.x = SPCDensity.X;
electrode.y = SPCDensity.Y;
electrode.z = SPCDensity.Z;
electrode.radius = SPCDensity.Size*22; % scale up just for visualization
electrode.color = hex2rgb(bands.color(SPCDensity.Color + 1));


fh{8} = figure('position',[200 200 400 400]);
cfg = [];
cfg.h_ax = gca;
cfg.view = [-20,15];
hold on

 %cfg.surface_facecolor = hex2rgb(bands.color(electrode.color + 1));
    %cfg.surface_facealpha = faceAlphas(target_i);

bml_plot3d_surface(cfg, vertices,faces);


hold on

cfg.h_ax = gca;
cfg.annotate = false;
bml_plot3d_points(cfg, electrode);







%
% 
% ROI_FrequencyBand_SPCdensity_perm = cell(6,numel(ROI_AtlasLabels));
% ROI_FrequencyBand_SPCdensity_perm_5pth = nan(6,numel(ROI_AtlasLabels));
% ROI_FrequencyBand_SPCdensity_perm_50pth = nan(6,numel(ROI_AtlasLabels));
% ROI_FrequencyBand_SPCdensity_perm_95pth = nan(6,numel(ROI_AtlasLabels));
% 
% 
% for fi = 2:6
%     tmp = ClusterEAtlasLabel(ClusterCentroidBand == fi);
%     % tmp = ClusterCentroidBand(strcmpi(ClusterEAtlasLabel,ROI_AtlasLabels{roi_i}));
%     for roi_i = 1 : numel(ROI_AtlasLabels)
%         for perm_i = 1 : nperms
%             tmp2 = ClusterEAtlasLabel(randperm(numel(ClusterEAtlasLabel), numel(tmp)));
% 
%             ROI_FrequencyBand_SPCdensity_perm{fi,roi_i}(perm_i) = 100*sum(strcmpi(tmp2,ROI_AtlasLabels{roi_i}))/ROI_nPairs(roi_i);
%         end
%         ROI_FrequencyBand_SPCdensity_perm_5pth(fi,roi_i) = prctile(ROI_FrequencyBand_SPCdensity_perm{fi,roi_i},5);
%         ROI_FrequencyBand_SPCdensity_perm_50pth(fi,roi_i) = prctile(ROI_FrequencyBand_SPCdensity_perm{fi,roi_i},50);
%         ROI_FrequencyBand_SPCdensity_perm_95pth(fi,roi_i) = prctile(ROI_FrequencyBand_SPCdensity_perm{fi,roi_i},95);
%     end
% end
% 
% fh{8} = figure('position',[200 200 1600 400]);
% tiledlayout(1,5)
% for fi = 2:6
%     nexttile
%     bar(1:numel(ROI_AtlasLabels),ROI_FrequencyBand_density(fi,:))
%     xticks(1:numel(ROI_AtlasLabels))
%     hold on
%      plot(1:numel(ROI_AtlasLabels),ROI_FrequencyBand_SPCdensity_perm_5pth(fi,:),'o')
%      plot(1:numel(ROI_AtlasLabels),ROI_FrequencyBand_SPCdensity_perm_50pth(fi,:),'o')
%    plot(1:numel(ROI_AtlasLabels),ROI_FrequencyBand_SPCdensity_perm_95pth(fi,:),'o')
% 
%     xticklabels(ROI_AtlasLabels)
%     ylabel('t-SPC density')
%     box off
%     ylim([0 15.2])
%     
% end
