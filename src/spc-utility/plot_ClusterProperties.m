function fh = plot_ClusterProperties(cfg, Clusters)
EvtTimes = cfg.EvtTimes;
Properties = cfg.Properties;
nProperties = numel(Properties);
bands = cfg.bands;
T = cfg.T;
fh = cell(1,nProperties);
% get ClusterCentroid
%ClusterCentroid = get_SPCClustersProperties(cfg.VarF,Clusters);
EvtTypes = cfg.cfg.plot.EvtTypes;
Task_labels = cfg.Task_labels;

for prop_i = 1 : nProperties
    Property = Properties{prop_i};
    switch Property
        case 'Frequency_1DHist'
            cfg.VarFields = {'ClusterCentroid'};
            ClusterCentroid = get_SPCClustersProperties(cfg,Clusters);

            fh{prop_i} = figure('Position',[400 200 600 400]);
            raincloud_plot(ClusterCentroid(:,2),"box_on",0,"color",[.8 .8 .8],"alpha",.3);
            hold on
            %xline(median(ClusterCentroid(:,2)),'r--',sprintf(['   %1.2f ' char(177)  ' %1.2f Hz'],median(ClusterCentroid(:,2)),iqr(ClusterCentroid(:,2))),'LabelOrientation','horizontal','linewidth',1.5)
            xline( median(ClusterCentroid(:,2)),'r--',sprintf('   %1.2f Hz ', median(ClusterCentroid(:,2))),'LabelOrientation','horizontal','linewidth',1.5)
            box off
            grid off
            xlabel(" Clusters f_{c} [Hz]" )
            ylabel(" PDF [a.u.] ")
            xticks([ 10 100])
            xlim([4 150])

            for i=1:height(bands)
                fill([bands.fstarts(i),bands.fstarts(i),bands.fends(i),bands.fends(i)],...
                    [-0.03,-0.025,-0.025,-0.03] - 0.01 , ...
                    hex2rgb(bands.color(i)),'EdgeColor','black','Marker','none','FaceAlpha',.85);
                text(bands.fmid(i)-1,-0.0275 - 0.01 ,bands.symbol{i});
            end
            set(gca,'xScale','log')


        case 'Duration_1DHist'
            % get ClusterDur
            cfg.VarFields = {'ClusterDur'};
            ClusterDur = get_SPCClustersProperties(cfg,Clusters);
            fh{prop_i} = figure('Position',[400 200 600 400]);
            raincloud_plot(ClusterDur,"box_on",0,"color",[.8 .8 .8],"alpha",.3)
            hold on
            xline(median(ClusterDur),'r--',sprintf(['   %1.3f ' char(177)  ' %1.3f s '],median(ClusterDur),iqr(ClusterDur)),'LabelOrientation','horizontal','linewidth',1.5)
            box off
            grid off
            xlabel(" Clusters Duration [s]" )
            ylabel(" PDF [a.u.] ")
            yticks([0 2 4])
            xlim([0 1])

        case 'TaskModulation_1DArea'
            cfg.VarFields = {'clust_prob','ClusterCentroidBand','Sign_Taskid'};
            [cluster_prob,ClusterCentroidBand,Sign_Taskid]  = get_SPCClustersProperties(cfg,Clusters);
            cluster_distr_bands = nan(6,numel(T));
            nClusters_all = size(cluster_prob,1);
            for ti = 1 : numel(T)
                tmp = ClusterCentroidBand(cluster_prob(:,ti)>0);
                for fi = 1:6
                    cluster_distr_bands(fi,ti) = sum(tmp ==  fi)/nClusters_all;
                end
            end
            fh{prop_i} =figure("renderer","painters","position",[400 200 600 400]);
            hold on
            area(T ,flipud(cluster_distr_bands)');
            colororder(flipud(hex2rgb(bands.color)))


            %legend(aa,bands.name,'location','best')
            set(gca, 'XTick', EvtTimes([2 3 4 5]))
            xLabeling = EvtTypes(2:end-1);
            set(gca, 'XTickLabel', xLabeling)
            xtickangle(30)
            ylabel( " t-SPC occurrence [a.u.] ")
            % distribution sample-by-sample frequencies
            xlim([-2.5 2])


        case 'TaskModulation_1DAggregationLine'

            cfg.VarFields = {'clust_prob','ClusterCentroidBand'};
            [cluster_prob,ClusterCentroidBand]  = get_SPCClustersProperties(cfg,Clusters);
            nClusters_all = size(cluster_prob,1);

            cluster_distr_bands = nan(6,numel(T));
            nClusters_all = size(cluster_prob,1);
            for ti = 1 : numel(T)
                tmp = ClusterCentroidBand(cluster_prob(:,ti)>0);
                for fi = 1:6
                    cluster_distr_bands(fi,ti) = sum(tmp ==  fi)/nClusters_all; % hardcode 2943
                end
            end

            cfg_cagg  = [];
            cfg_cagg.nperms = 500;
            cfg_cagg.T = T;
            cluster_aggr_bands = nan(6,numel(T));
            cluster_probperm_bands = nan(6, numel(T), cfg_cagg.nperms);

%             if isfile('permClusters.mat')
%                 load('permClusters.mat')
%                 for fi = 2:6
%                     cluster_probperm_bands(fi,:,:) = cluster_probperm_bands(fi,:,:)*sum(ClusterCentroidBand == fi)/nClusters_all;
%                 end
%             else
                for fi = 2:6
                    [cluster_aggr_bands(fi,:),cluster_probperm_bands(fi,:,:)] = compute_clusteraggregation(cfg_cagg, cluster_prob(ClusterCentroidBand == fi,:));
                    cluster_probperm_bands(fi,:,:) = cluster_probperm_bands(fi,:,:)*sum(ClusterCentroidBand == fi)/nClusters_all;
                end

                [cluster_aggr, cluster_prob_perm] = compute_clusteraggregation(cfg_cagg, cluster_prob);
%             end

            cluster_distr_bands_task = nan(6,5);
            for evt = 1:5
                cluster_distr_bands_task(:,evt) = mean(cluster_distr_bands(:,T <= EvtTimes(:,evt+1) & T >= EvtTimes(:,evt)),2,'omitnan');
            end

            cluster_distr_task = sum(cluster_distr_bands_task);

            % plot occurrence and permutation
            cluster_prob_perm_task = nan(cfg_cagg.nperms,5);
            cluster_prob_perm_bandtask = nan(cfg_cagg.nperms,6,5);
            for evt = 1:5
                cluster_prob_perm_task(:,evt) = mean(cluster_prob_perm(T <= EvtTimes(evt+1) & T >= EvtTimes(evt),:),1,'omitnan');
                cluster_prob_perm_bandtask(:,:,evt) = squeeze(mean(cluster_probperm_bands(:,T <= EvtTimes(evt+1) & T >= EvtTimes(evt),:),2,'omitnan'))';
            end

            % compute index
            index_cluster_aggr = mean(cluster_prob) - median(cluster_prob_perm,2,'omitnan')';

            fh{prop_i} = figure('position',[200 200 500 1200]);
            tiledlayout(6,1)
            nexttile
            plot(T,mean(cluster_prob))
            hold on
            plotShaded(T,[prctile(cluster_prob_perm,95,2)';  median(cluster_prob_perm,2)'; prctile(cluster_prob_perm,5,2)'], [0 0 0],'-',.2)
            xlim([T(1) T(end)])
            box off
            ylabel(' Total ')

            for fi = 2:6
                nexttile
                plot(T,100*cluster_distr_bands(fi,:))
                hold on
                plotShaded(T,100*[squeeze(prctile(cluster_probperm_bands(fi,:,:),5,3)); ...
                    squeeze(prctile(cluster_probperm_bands(fi,:,:),50,3)); ...
                    squeeze(prctile(cluster_probperm_bands(fi,:,:),95,3))], [0 0 0], '-',.2)
                xlim([T(1) T(end)])
%                 if  fi == 3 || fi == 5 || fi == 6
%                     ylim([0 0.033])
%                 end
                %ylim([0 0.08])
                box off
                ylabel(bands.symbol(fi))
            end


            
      
        case 'TaskModulation_1DPhaseAggregation'
            cfg.VarFields = {'clust_prob','ClusterPhase','ClusterOnOff','ClusterCentroidBand'};
            [cluster_prob, ClusterPhase, ClusterOnOff, ClusterCentroidBand] = get_SPCClustersProperties(cfg,Clusters);
            nClusters_all = size(cluster_prob,1);

            cluster_phase = nan(nClusters_all,numel(T));
            for ii = 1 : nClusters_all
                cluster_phase(ii,T<= (ClusterOnOff(ii,2)) & T >= (ClusterOnOff(ii,1))) = ClusterPhase(ii,1);
            end

            cluster_distr_bands = nan(6,numel(T));
            for ti = 1 : numel(T)
                tmp = ClusterCentroidBand(cluster_prob(:,ti)>0);
                for fi = 1:6
                    cluster_distr_bands(fi,ti) = sum(tmp ==  fi)/nClusters_all;
                end
            end



            winmov_size = round(0.005/0.005);%round(0.2/clust_res);
            winmov_step = round(0.005/0.005);

            winmov = bsxfun(@plus, (0:winmov_step:(numel(T)-winmov_size))', 1:winmov_size);
            nwinmov = size(winmov,1);
            ClusterPhaseinTrial = cell(6,nwinmov);
            ClusterPhaseinTrial_confmean = nan(3,6,nwinmov);
            ClusterPhaseinTrial_mvlstat = nan(6,nwinmov);
            ClusterPhaseinTrial_pvalstat = nan(6,nwinmov);
             ClusterPhaseinTrial_stdstat = nan(6,nwinmov);
           
            %ClusterPhaseinTrial_circstd = 
            for wini = 1 : nwinmov
                for kk = 1 : winmov_size
                    tmp = cluster_prob(:,winmov(wini,kk))>0;
                    for fi = 2:6
                        ClusterPhaseinTrial{fi,wini} = cluster_phase(ClusterCentroidBand == fi & tmp',winmov(wini,kk));
                        ClusterPhaseinTrial{fi,wini} = ClusterPhaseinTrial{fi,wini}(~isnan(ClusterPhaseinTrial{fi,wini} ));
                    end
                end
            end



            for wini = 1 : nwinmov
                for fi = 2:6
                    if numel(ClusterPhaseinTrial{fi,wini} ) > 3
                        [ClusterPhaseinTrial_confmean(1,fi,wini),ClusterPhaseinTrial_confmean(2,fi,wini),ClusterPhaseinTrial_confmean(3,fi,wini)] = circ_mean(ClusterPhaseinTrial{fi,wini});
                        ClusterPhaseinTrial_pvalstat(fi,wini) = circ_otest(ClusterPhaseinTrial{fi,wini});

                        ClusterPhaseinTrial_mvlstat(fi,wini) = circ_r(ClusterPhaseinTrial{fi,wini});
                        ClusterPhaseinTrial_stdstat(fi,wini) = circ_std(ClusterPhaseinTrial{fi,wini});
                    else
                        ClusterPhaseinTrial_confmean(:,fi,wini) = nan(1,3);
                        ClusterPhaseinTrial_pvalstat(fi,wini) = nan;
                        ClusterPhaseinTrial_mvlstat(fi,wini) = nan;
                        ClusterPhaseinTrial_stdstat(fi,wini) = nan;                    
                    end
                end
            end
            %
            Twin = T(winmov(:,ceil(winmov_size/2)));

            fh{prop_i} = figure("renderer","painters","position",[400 200 500 1200]);
            tiledlayout(5,1)
            for fi = 2:6
                nexttile

                % nonanplot = ~isnan(ClusterPhaseinTrial_confmean(1,fi,:));
                %     errorbar(Twin(nonanplot),squeeze(ClusterPhaseinTrial_confmean(1,fi,nonanplot)),squeeze(ClusterPhaseinTrial_confmean(1,fi,nonanplot)) - squeeze(ClusterPhaseinTrial_confmean(3,fi,nonanplot)), ...
                %         squeeze(ClusterPhaseinTrial_confmean(1,fi,nonanplot)) - squeeze(ClusterPhaseinTrial_confmean(2,fi,nonanplot)), ...
                %         'o','color',[0 0 0 .4],'MarkerSize',.1,'LineWidth',0.01,...
                %     'MarkerEdgeColor',[0 0 0 ],'MarkerFaceColor','k')
                scatter(Twin(nonanplot),squeeze(ClusterPhaseinTrial_confmean(1,fi,nonanplot)),8, ...
                    'ok','filled')
                  
                hold on
%                 errorbar(Twin(nonanplot), squeeze(ClusterPhaseinTrial_confmean(1,fi,nonanplot)), ... squeeze(ClusterPhaseinTrial_confmean(2,fi,nonanplot)),...
%                     ClusterPhaseinTrial_stdstat(fi,nonanplot) )
                %     plot(Twin,squeeze(ClusterPhaseinTrial_confmean(2,fi,:)))
                %     plot(Twin,squeeze(ClusterPhaseinTrial_confmean(3,fi,:)))
                xlim([-2.5 2])


                for evt = [2 3 4 5]'
                    plot(EvtTimes(evt)' *[1,1], [-pi pi], '--', 'Color', 'r', 'LineWidth', 1)
                end
                set(gca, 'XTick', EvtTimes([2 3 4 5]))
                xLabeling = EvtTypes(2:end-1);
                %xLabeling{2} = '';
                set(gca, 'XTickLabel', xLabeling)
                xtickangle(30)
                yticks([-pi -pi/2 0 pi/2 pi])
                set(gca,'YTickLabel',{'-\pi','-\pi/2','0','\pi/2','\pi'})
                h = fdr_bh(ClusterPhaseinTrial_pvalstat(fi,:))>0;
                h = ClusterPhaseinTrial_pvalstat(fi,:) <= .05;
                hh = bwconncomp(h);
                for xx = 1  : hh.NumObjects
                    patch([Twin(hh.PixelIdxList{xx}(1)) Twin(hh.PixelIdxList{xx}(end)) Twin(hh.PixelIdxList{xx}(end)) Twin(hh.PixelIdxList{xx}(1))],[4.7 4.7 5.1 5.1],'k','facealpha',1)
                end
                %patch(Twin(h == 1) , 5*ones(1,sum(h == 1)),'k*')
                box off
                ylim([-5/4*pi 5.1])
                ylabel("\phi_{pref}")
                title(bands.symbol(fi))
                yyaxis right
                plot(T, cluster_distr_bands(fi,:)*nClusters_all,'linewidth',.7,'color',hex2rgb(bands.color(fi)))
                ylabel('Number Clusters')
                set(gca,'YColor',hex2rgb(bands.color(fi)))
                %ylim([ 0 30])
            end


            %             tiledlayout(2,1)
            %             nexttile
            %             polarhistogram(ClusterPhase,18,'FaceColor',[.8 .8 .8],'normalization','pdf')
            %
            %             thetaticks([0 90 180 270])
            %             thetaticklabels({'0','\pi/2','\pm \pi','-\pi/2'})
            %
            %             set(gca,'fontsize',18)
            %
            %             nComponents = 1;
            %             fittedVmm = fitmvmdist(ClusterPhase(:,1), nComponents, ...
            %                 'MaxIter', 250); % Set maximum number of EM iterations to 250
            %             hold on
            %             % Plot initial and fitted distributions.
            %             likelihoodsFitted = fittedVmm.pdf(linspace(-pi,pi,1000)');
            %             polarplot(linspace(-pi,pi,1000)', likelihoodsFitted);
            %             [p,z] = circ_otest(ClusterPhase(:,1));
            %             title(sprintf("Hodges-Ajne test: p = %1.2f ", p))
            %             %axis([-pi, pi, 0, 1]);
            %             % build clusters phase
            %
            %             cluster_phase = nan(nClusters_all,numel(T));
            %             cluster_ES = nan(nClusters_all,numel(T));
        case 'TaskModulation_LagDelay'
            cfg.VarFields = {'ClusterPhase','ClusterCentroid'};
            [ClusterPhase, ClusterCentroid] = get_SPCClustersProperties(cfg,Clusters);

            phase = ClusterPhase(:,1);
            freqs = ClusterCentroid(:,2);
            nbtsp = 500;
            FreqRange = cfg.freq_range;
            idxLow = (freqs <= FreqRange(2) & freqs >= FreqRange(1));
            % between 4-30
            [rLow,delayLow, pvalLow, fLow] = compute_delay_withperm(freqs(idxLow) ,phase(idxLow),nbtsp);
            set(fLow,'visible','on')


            fh{prop_i} = fLow;

        case 'TaskModulation_Phasepolar'
            cfg.VarFields = {'ClusterPhase','ClusterCentroid','PPCz'};
            [ClusterPhase, ClusterCentroid,PPCz] = get_SPCClustersProperties(cfg,Clusters);
            FreqRange = cfg.freq_range;
            idxLow = (freqs <= FreqRange(2) & freqs >= FreqRange(1));
            phase = ClusterPhase(:,1);
            freqs = ClusterCentroid(:,2);
            fh{prop_i}  = polarplot_phase_freq(phase(idxLow), PPCz(idxLow), freqs(idxLow), 55, hex2rgb([bands.color(2:3); '#ffef5c']));
      
        case 'TaskModulation_Lagdelay_overtime'

          % take (theta+alpha)/beta ratio sPC
          cfg.VarFields = {'ClusterPhase','ClusterCentroid','clust_prob','ClusterCentroidBand'};
            [ClusterPhase, ClusterCentroid,cluster_prob,ClusterCentroidBand] = get_SPCClustersProperties(cfg,Clusters);
            nClusters_all = size(cluster_prob,1);
            cluster_distr_bands = nan(6,numel(T));
            nClusters_all = size(cluster_prob,1);
            for ti = 1 : numel(T)
                tmp = ClusterCentroidBand(cluster_prob(:,ti)>0);
                for fi = [2 3 4]
                    cluster_distr_bands(fi,ti) = sum(tmp ==  fi)/nClusters_all; % hardcode 2943
                end
            end
            
            %cluster_distr_ratio = sum(cluster_distr_bands(2:3,:))./cluster_distr_bands(4,:);
            cluster_distr_ratio = (sum(cluster_distr_bands(2:3,:)) - cluster_distr_bands(4,:))./sum(cluster_distr_bands(2:4,:));

                [rLow_time, delayLow_time, pvalLow_time] = deal(nan(1,numel(T)));

            FreqRange = cfg.freq_range;
            idxLow = (freqs <= FreqRange(2) & freqs >= FreqRange(1));
            phase = ClusterPhase(:,1);
            freqs = ClusterCentroid(:,2);
                for ti = 1 : numel(T)
                    freqs_ = freqs(cluster_prob(:,ti) == 1);
                    phase_ = phase(cluster_prob(:,ti) == 1);
                    freqs__ = freqs(cluster_prob(:,ti) == 1 & idxLow);
                    phase__ = phase(cluster_prob(:,ti) == 1 & idxLow);
                    try
                        [rLow_time(ti),delayLow_time(ti), pvalLow_time(ti)] = compute_delay_withperm(freqs__ ,phase__, nbtsp);
                    catch
                        warning('error time phase computation: leave NaN value')
                    end
                end
                rLow_time(pvalLow_time > .05) = nan; % nan

                
                fh{prop_i} = figure('renderer','painters','Position',[200 200 900 800]);
                nexttile
                scatter(T(pvalLow_time < .05),delayLow_time(pvalLow_time < .05)*1E+3,15,'k','filled')
                hold on
                for evt = [2 3 4 5]'
                    plot(EvtTimes(evt)' *[1,1], [min(delayLow_time) max(delayLow_time)]*1E+3, '--', 'Color', 'w', 'LineWidth', 1)
                end
                set(gca, 'XTick', EvtTimes([2 3 4 5]))
                yline(0,'--')
                xLabeling = EvtTypes(2:end-1);
                set(gca, 'XTickLabel', xLabeling)
                xtickangle(30)
                ylabel('Time delay [ms]')
                box off
                ylim([-120 120])

                nexttile
                scatter(T,rLow_time,15,'filled')
                hold on

                for evt = [2 3 4 5]'
                    plot(EvtTimes(evt)' *[1,1], [min(delayLow_time) max(delayLow_time)]*1E+3, '--', 'Color', 'w', 'LineWidth', 1)
                end
                set(gca, 'XTick', EvtTimes([2 3 4 5]))
                yline(0,'--')
                xLabeling = EvtTypes(2:end-1);
                set(gca, 'XTickLabel', xLabeling)
                xtickangle(30)
                ylabel('R')
                box off
                ylim([-1 1])
                yyaxis right
                plot(T,cluster_distr_ratio)
          


               
% % over time
% [rAll_time, delayAll_time, pvalAll_time] = deal(nan(1,numel(T)));
% 

% 
% rAll_time(pvalAll_time > .05) = nan;
% rLow_time(pvalLow_time > .05) = nan;
% %
% figure('renderer','painters','Position',[200 200 900 300])
% nexttile
% scatter(T(pvalAll_time < .05),delayAll_time(pvalAll_time < .05),15,'k','filled')
% 
% nexttile
% plot(T,sum(cluster_prob,1),'linewidth',1.5)
% nexttile
% plot(T,-log(pvalAll_time),'linewidth',1.5)
% hold on
% yline(-log(.05))
% nexttile
% plot(T,rAll_time,'linewidth',1.5)


            
    end

end






%% temporary material




%%
% figure('renderer','painters','Position',[200 200 900 400])
% nexttile
% scatter(T(pvalLow_time < .05),delayLow_time(pvalLow_time < .05)*1E+3,15,'k','filled')
% hold on
% for evt = [2 3 4 5]'
%     plot(meanEvts(evt)' *[1,1], [min(delayLow_time) max(delayLow_time)]*1E+3, '--', 'Color', 'w', 'LineWidth', 1)
% end
% set(gca, 'XTick', meanEvts([2 3 4 5]))
% xLabeling = cfg.plot.EvtTypes(2:end-1);
% set(gca, 'XTickLabel', xLabeling)
% xtickangle(30)
% ylabel('Time delay [ms]')
% box off
% 
% nexttile
% plot(T,sum(cluster_prob,1),'linewidth',1.5)
% box off
% hold on
% for evt = [2 3 4 5]'
%     plot(meanEvts(evt)' *[1,1], [0 max(sum(cluster_prob,1))], '--', 'Color', 'w', 'LineWidth', 1)
% end
% set(gca, 'XTick', meanEvts([2 3 4 5]))
% xLabeling = cfg.plot.EvtTypes(2:end-1);
% set(gca, 'XTickLabel', xLabeling)
% xtickangle(30)
% nexttile
% plot(T,-log(pvalLow_time),'linewidth',1.5)
% box off
% hold on
% for evt = [2 3 4 5]'
%     plot(meanEvts(evt)' *[1,1], [0 max(-log(pvalLow_time))], '--', 'Color', 'w', 'LineWidth', 1)
% end
% set(gca, 'XTick', meanEvts([2 3 4 5]))
% xLabeling = cfg.plot.EvtTypes(2:end-1);
% set(gca, 'XTickLabel', xLabeling)
% xtickangle(30)
% yline(-log(.05))
% 
% nexttile
% plot(T,rLow_time,'linewidth',1.5)
% hold on
% box off
% for evt = [2 3 4 5]'
%     plot(meanEvts(evt)' *[1,1], [0 max(rLow_time)], '--', 'Color', 'b', 'LineWidth', 1)
% end
% set(gca, 'XTick', meanEvts([2 3 4 5]))
% xLabeling = cfg.plot.EvtTypes(2:end-1);
% set(gca, 'XTickLabel', xLabeling)
% xtickangle(30)
% 
% %%
% figure('renderer','painters','Position',[200 200 900 400])
% nexttile
% scatter(T(pvalLow_time < .05),delayLow_time(pvalLow_time < .05)*1E+3,15,'k','filled')
% hold on
% for evt = [2 3 4 5]'
%     plot(meanEvts(evt)' *[1,1], [min(delayLow_time) max(delayLow_time)]*1E+3, '--', 'Color', 'w', 'LineWidth', 1)
% end
% set(gca, 'XTick', meanEvts([2 3 4 5]))
% yline(0,'--')
% xLabeling = cfg.plot.EvtTypes(2:end-1);
% set(gca, 'XTickLabel', xLabeling)
% xtickangle(30)
% ylabel('Time delay [ms]')
% box off
% ylim([-120 120])
% yyaxis right
% scatter(T,rLow_time,15,'filled')
% ylabel(' R ')
% 
% 
% figname = fullfile(PATHS.saveFigures,'delayLow_overtime');
% saveas(gcf,[figname, '.fig'])
% saveas(gcf,[figname, '.pdf'])
% saveas(gcf,[figname, '.svg'])
% saveas(gcf,[figname, '.png'])
% mean(pvalLow_time <=  .05 & rLow_time > 0)
% mean(pvalLow_time <=  .05 & rLow_time < 0)
% 
% %% try polar plot
% 
% figure('Position',[300 300 400 400])
% fh = polarplot_phase_freq(phase(idxLow), PPCz(idxLow), freqs(idxLow), 55);
% figname = fullfile(PATHS.saveFigures,'delayLow_relationship_polar');
% saveas(fh,[figname, '.fig'])
% saveas(fh,[figname, '.pdf'])
% saveas(fh,[figname, '.svg'])
% saveas(fh,[figname, '.png'])