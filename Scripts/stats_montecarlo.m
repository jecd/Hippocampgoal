function [ mc_stats ] = stats_montecarlo( cfg, data1,data2 )
%stats_montecarlo - nonparametric statistical test by calculating
% Monte-Carlo estimates of the significance probabilities and/or critical values
% from the permutation distribution
%INPUTS
% cfg - conficuration struct containing various user specified options 
% data1...N - data to statistically compare
%OUTPUTS
% stats = structure containing stats for comparison
% JCD adapted from KK 


D = data1 - data2; % get difference mat

% faster than matlab defined functions 
real_t = (squeeze(mean(D, 1))) ./ (squeeze((std(D, 0, 1))/sqrt(size(D, 1)))); % SUBJECT is in dimension1 in this data


% preallocate
perm_tmap = nan(size(D, 2),  size(D, 3), cfg.nperm); % nTR * nTR * perm 
max_clust_t_neg  = zeros(cfg.nperm,1); % storage  
max_clust_t_pos  = zeros(cfg.nperm,1); % storage  
% get permed t-maps
for iPerm = 1: 1: cfg.nperm
    % compute t-map of null hypothesis
    rng('shuffle');
    tmp_perm = (double(rand(size(D, 1), 1) > 0.5)-0.5)*2; % this shuffels data at the condition level
    perm_mult = repmat(reshape(tmp_perm, [size(D,1), 1, 1]),  [1 size(D, 2), size(D, 3)]);% generate the randomly shuffled matrix
    tmp_shuf = D.*perm_mult; % applt shuffle 
    tnum = squeeze(mean(tmp_shuf, 1)); % get mean diff
    tdenom = squeeze(std(tmp_shuf, 0, 1)) /sqrt(size(D, 1)-1); % sem calc
    tmap   = tnum./tdenom;

    perm_tmap (:, :, iPerm) = tmap; % save 
    tmap(abs(tmap)<tinv(1-(cfg.alpha/2), size(D, 1)-1))=0; % 2tailed

    
    % grab clusters & t sums for each cluster
    tmap = tril(tmap); 
    clustinfo = bwconncomp(tmap); % find clusters
    tmp_mass = cell2mat([0 cellfun(@(x) sum(tmap(x)), clustinfo.PixelIdxList, 'UniformOutput', false)]); % sum clsuters 
    tmp_neg = min(tmp_mass(tmp_mass < 0)); % get max 
    tmp_pos = max(tmp_mass(tmp_mass >= 0)); % get max
    % ask kamin
    if isempty(tmp_neg) % for when there are no negative clusters
        max_clust_t_neg(iPerm) = 0;
    else
        max_clust_t_neg(iPerm) = tmp_neg;
    end
    if isempty(tmp_pos) % for when there are no positive clusters
        max_clust_t_pos(iPerm) = 0;
    else
        max_clust_t_pos(iPerm) = tmp_pos;
    end    
end 

% get z-score
zmap = (real_t-squeeze(mean(perm_tmap,3)))./squeeze(std(perm_tmap, 0, 3));

% apply cluster-level corrected threshold
zmapthresh = zmap;
zmapthresh(abs(zmapthresh)<norminv(1-cfg.alpha/2))=0; 

% find islands and remove those smaller than cluster size threshold
zmapthresh = tril(zmapthresh);  
clustinfo = bwconncomp(zmapthresh); % find clusters
clust_ts = cell2mat([cellfun(@(x) sum(zmapthresh(x)), clustinfo.PixelIdxList, 'UniformOutput', false)]); % sum clusters for mass
% threshhold seperately for pos and neg
% pos 
pos_clust_t_threshold = prctile(max_clust_t_pos, (1-cfg.clust_alpha)*100); % here I used mcc_cluster_pval = .05; 
% neg 
neg_clust_t_threshold = prctile(max_clust_t_neg, (cfg.clust_alpha)*100); % here I used mcc_cluster_pval = .05; 

% identify clusters to remove
% whichclusters2remove = find(clust_ts > neg_clust_t_threshold);

% identify sig clusters
sig_neg = find(clust_ts < neg_clust_t_threshold);
sig_pos = find(clust_ts > pos_clust_t_threshold);

% make an output struct 
mc_stats.pos_dist = max_clust_t_pos;
mc_stats.neg_dist = max_clust_t_neg;
mc_stats.tmap = real_t;
mc_stats.zmap = zmap;
mc_stats.zmapthresh = zmapthresh;

% assign neg clusters to struct
if ~isempty(sig_neg)
    sigmask = tril(zeros(size(D,2),size(D,3)));
    for iclust = 1:length(sig_neg)
        % make mask
        sigmask(clustinfo.PixelIdxList{sig_neg(iclust)}) = iclust; % linear index is faster
        % save things out
        mc_stats.neg_clust(iclust).pixels = clustinfo.PixelIdxList{sig_neg(iclust)}; % pixels
        mc_stats.neg_clust(iclust).clustermass = clust_ts(sig_neg(iclust)); % clustermass
        mc_stats.neg_clust(iclust).pval = length(find(clust_ts(sig_neg(iclust)) > max_clust_t_neg))/length(max_clust_t_neg + 1); % pval (number of clusters bigger than my cluster)/number of perms + 1
    end
    mc_stats.neg_sigmask = sigmask; % for visualizing 
else
    
end

% assign pos clusters to struct
if ~isempty(sig_pos)
    sigmask = tril(zeros(size(D,2),size(D,3)));
    for iclust = 1:length(sig_pos)
        % make mask
        sigmask(clustinfo.PixelIdxList{sig_pos(iclust)}) = iclust; % linear index is faster
        % save things out
        mc_stats.pos_clust(iclust).pixels = clustinfo.PixelIdxList{sig_pos(iclust)}; % pixels
        mc_stats.pos_clust(iclust).clustermass = clust_ts(sig_pos(iclust)); % clustermass
        mc_stats.pos_clust(iclust).pval = length(find(clust_ts(sig_pos(iclust)) < max_clust_t_pos))/length(max_clust_t_pos + 1); % pval (number of clusters bigger than my cluster)/number of perms + 1        
    end
    mc_stats.pos_sigmask = sigmask; % for visualizing 
else
    
end


end

