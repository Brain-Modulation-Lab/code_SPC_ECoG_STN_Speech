function [ROI_AtlasLabels,E_ROIAtlasLabel_pairs_n]  = getROIs(ECoG_atlas, ECoG_subject, min_subj,min_pairs,flag_sort)


AtlasLabels = unique(ECoG_atlas);

E_AtlasLabel_pairs_n = cellfun(@(x) sum(strcmpi(ECoG_atlas,x)),AtlasLabels);
E_AtlasLabel_subj_n = cellfun(@(x) numel(unique(ECoG_subject(strcmpi(ECoG_atlas,x)))),AtlasLabels);

ROI_AtlasLabels = AtlasLabels(E_AtlasLabel_pairs_n >= min_pairs & E_AtlasLabel_subj_n >= min_subj);
E_ROIAtlasLabel_pairs_n = E_AtlasLabel_pairs_n(E_AtlasLabel_pairs_n >= min_pairs & E_AtlasLabel_subj_n >= min_subj);

if ~exist("flag_sort",'var')
    flag_sort = false;
end

if flag_sort
    [~, roi_order] = sort(E_ROIAtlasLabel_pairs_n,'descend');
    ROI_AtlasLabels = ROI_AtlasLabels(roi_order);
    E_ROIAtlasLabel_pairs_n = E_ROIAtlasLabel_pairs_n(roi_order);
end