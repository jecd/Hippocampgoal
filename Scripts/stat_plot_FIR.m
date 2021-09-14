function [ main_fig ] = stat_plot_FIR( data1,data2,stats )
%stat_plot_FIR - this function plots stats for timepoint*timepoint fir correlation analysis
%INPUTS
% data1 - non subtracted data
% data2 - non subtracted data
% stats - some sort of stats structure. Must have field tmap or zmap and sigmask
%OUTPUTS 
% main_fig - a figure with format specified below

% make the square figure
main_fig=figure('Position', [0, 300, 1000, 1000],'visible','on');
subplot(2,2,1);d1_h=imagesc(squeeze(mean(data1,1)));d1_ax=gca;title(deunderscore(stats.cName1));colorbar;
subplot(2,2,2);d2_h=imagesc(squeeze(mean(data2,1)));d2_ax=gca;title(deunderscore(stats.cName2));colorbar;
% make the clims the same
d2_ax.CLim = d1_ax.CLim;

% plot the stats
if isfield(stats,'zmap') % try zmap
    subplot(2,2,3);imagesc(stats.zmap);colorbar;title('analytic zmap')
elseif isfield(stats,'tmap') % try tmap
end

if isfield(stats,'zmapthresh')
    subplot(2,2,4);imagesc(stats.zmapthresh>0 | stats.zmapthresh<0 );ax = gca;ax.CLimMode = 'manual';hold on;title('sig mask')
    
    if isfield(stats, 'neg_sigmask') % plot surviving clusters % might need to edit this so will plot contours better
        [lines,hlines]=contour(stats.neg_sigmask,1);ax_lines = gca;hlines.LineColor = 'cyan';hlines.LineWidth=2;
    end
    if isfield(stats, 'pos_sigmask')
        [lines,hlines]=contour(stats.pos_sigmask,1);ax_lines = gca;hlines.LineColor = 'cyan';hlines.LineWidth=2;
    end
elseif isfield(stats,'tmapthresh')
end

end

