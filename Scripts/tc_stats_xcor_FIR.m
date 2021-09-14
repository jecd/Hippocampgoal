function [ ] = tc_stats_xcor_FIR( data, contrasts, outdir )
%tc_stats_xcor_FIR - this function does paired t-tests between timepoint*timepoint corellation matrices 
%INPUTS - 
% data - data structure organized ROIs.mean_data.conditions
% contrasts - cell of names of conditions to compare. Will always do 1 vs 2. Supports multiple conditions
% OUTPUTS - will save a ROI specific directory structure within outdir

ROI_names = fieldnames(data);
contrast_names = fieldnames(data.(ROI_names{1}).mean_data);

% make plots dir
if ~exist(outdir, 'dir')
    mkdir(outdir)
end

for iROI = 1:length(ROI_names)
    cur_roi = ROI_names{iROI}; % name for index
    
    % make plots dir if it doesnt exist
    plot_dir = [outdir, cur_roi, '/'];
    if ~exist(plot_dir, 'dir')
        mkdir(plot_dir)
    end
    
    for icontrast = 1:length(contrasts)
        % names for index
        cName1 = contrasts{icontrast}{1};
        cName2 = contrasts{icontrast}{2};
        
        % do the index
        data1 = data.(cur_roi).mean_data.(cName1).indiv;
        data2 = data.(cur_roi).mean_data.(cName2).indiv;
        
        % run the contrast
        disp(['running ', cName1, ' vs ' cName2 ])
        disp('permuting ')
        stats = [];
        cfg = [];
        cfg.nperm = 10000;
        cfg.stat= 'Tstat';
        cfg.alpha = 0.05; % this will be divided by two because two.tailed test.
        cfg.clust_alpha = 0.05; % this is for the cluster
        [stats]=stats_montecarlo(cfg,data1,data2);
        % add names in for plotting 
        stats.cName1 = cName1;
        stats.cName2 = cName2;
        stats.ROI = cur_roi;
        % add in data for exporting. Names are cName1=data1 and cName2=data2
        stats.data1 = data1;
        stats.data2 = data2;
                
        % visualize 
        disp('plotting...')
        [fig2save]=stat_plot_FIR(data1,data2,stats);

        % save 
        saveas(fig2save,[plot_dir, cur_roi, '_', cName1, '_vs_', cName2, '.eps'], 'epsc')
        % close
        close(fig2save)
        
        % save
        all_stats(icontrast).stats = stats;
                
    end % compare
    %% save 
    save([plot_dir,cur_roi, '_stats.mat'], 'all_stats')
    
    clear all_stats
end% roi

end

