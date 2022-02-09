close all; clear;

addpath("Colormaps/")
%% Plot dir
plot_dir = '/Volumes/Data/zoocon/Hippocampgoal/Figures/';
%% Set parameter
% gamma=.01;              % decay rate of successor representatin (SR) 
gamma=.3;              % decay rate of successor representatin (SR) 
lambda=gamma;      % decay rate of expected value (EV)

%% Create Zoo maze
Zoo_s=[1,2,3,4,6,7,3,8];
Zoo_t=[2,3,4,5,7,3,8,9];
Zoo=graph(Zoo_s, Zoo_t); % no directional graph
Zoo_adj=full(adjacency(Zoo)); % adjacency matrix

ZooStr=figure;
subplot(1,3,1); plot(Zoo, 'Layout','force'); title('Structure');
subplot(1,3,2); imagesc(Zoo_adj); title('Adjacency');

%% Get the Path
thepath=shortestpath(Zoo,1,5); % The original path 1->5
path_div=shortestpath(Zoo,1,9); % A diverging path 1->9
path_con=shortestpath(Zoo,6,5); % A converging path 6->5 

%% Compute SR
Mat=Zoo_adj;
zSR1=inv(eye(size(Mat))-gamma*Mat); % SR matrix
subplot(1,3,3); imagesc(zSR1); title(['SR gamma=',  num2str(gamma)]);
colormap(viridis)
set(gcf, 'position', [100 100 700 200]);
saveas(gcf, [plot_dir, 'Fig_S1_zoo_SR_gamma_', num2str(gamma),'.eps'],'epsc');
close

% SR at each of 5 points of navigation 
Cd_path_bi=zSR1(:,thepath); % SR when participants nagivating 1->5
Cd_div_bi=zSR1(:,path_div); % SR when participants nagivating 1->9
Cd_con_bi=zSR1(:,path_con); % SR when participants nagivating 6->5

%% Compute sim between SR
PS_same_H1 = corr(Cd_path_bi(:,1), Cd_path_bi(:,1));
PS_div_H1 = corr(Cd_path_bi(:,1), Cd_div_bi(:,1));
PS_con_H1 = corr(Cd_path_bi(:,1), Cd_con_bi(:,1));
PS_no_overlap_H1 = corr(Cd_div_bi(:,1), Cd_con_bi(:,1));

PS_Fig=figure; 
bar([PS_same_H1,PS_con_H1, PS_div_H1, PS_no_overlap_H1])
xticklabels({'Same Sequence', 'Converging','Diverging', 'No Overlap'})
ylabel('Similarity');
ax = gca;
set(ax,'TickDir','out','box','off')
ylim([-1,1.1])
title('Pattern Similarity (Pearsons)')
set(gcf, 'position', [100 100 700 200]);
saveas(gcf, [plot_dir, 'Fig_S1_zoo_SR_pattern_similarity_gamma_', num2str(gamma),'.eps'],'epsc');
close 