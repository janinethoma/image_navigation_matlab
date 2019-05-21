%% Demo for Mapping, Localization and Path Planning for Image-based Navigation using Visual Features and Map
%
% This demo does the following using the concepts introduced in our paper:
% 1) Find landmarks in a sub sequence of the Oxford Robotcar run from  2015-10-29 12:18:17
% 2) Match a short query sequence from 2014-11-18 13:20:12
% Both, reference and query, had to be shortened due to supplementary
% material data limis. 
%
% If you do not have mosek installed, you can have a look at the saved
% figures in the results folder instead.

%% CHANGE THIS TO MATCH YOUR MOSEK INSTALLATION
dpath = javaclasspath;
if isempty(dpath)
    path = '/path/to/your/mosek/installation';
    system(['export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:', path, '/mosek/8/tools/platform/linux64x86/bin/libmosek64.so.8.0']);
    addpath([path, '/mosek/8/toolbox/r2014a/']);
    javaaddpath([path, '/mosek/8/tools/platform/linux64x86/bin/mosekmatlab.jar']);
end


addpath('helper')

%% Load image locations and precalculated NetVLAD image feature distances

load('data.mat', 'ref', 'query', 'f_dists_ref_ref', 'f_dists_query_ref');

ref.x = ref.xy(1,:);
ref.y = ref.xy(2,:);

query.x = query.xy(1,:);
query.y = query.xy(2,:);

f1 = figure; hold on;
scatter(ref.x, ref.y,'k','.');
scatter(query.x, query.y,[],query.t,'*');
legend({'Reference', 'Query'});
title('Reference and query sequence');
xlabel('UTM Easting [m]')
ylabel('UTM Northing [m]')
hold off;
saveas(f1, fullfile('results', 'original_sequences.eps'),'epsc')


%% Generate topology
[edges, geo_dists] = get_edges(ref.xy,10,2);
source_idx = 1;
sink_idx = length(ref.x);

f2 = figure; hold on;
for i=1:50:length(geo_dists)
    P1 = ref.xy(:,edges(1,i));
    P2 = ref.xy(:,edges(2,i));
    plot([P2(1);P1(1)],[P2(2);P1(2)], 'HandleVisibility','Off');
    hold on;
end

scatter(ref.x(source_idx), ref.y(source_idx),'r','o');
scatter(ref.x(sink_idx), ref.y(sink_idx),'g','o');
legend({'Source','Sink'})
title('Topology')
xlabel('UTM Easting [m]')
ylabel('UTM Northing [m]')
hold off;
saveas(f2, fullfile('results', 'topology.eps'),'epsc')


%% Sample landmarks with flow
num_landmarks = 100;

% Get feature distance for each edge in reference topology
feat_dists = f_dists_ref_ref(sub2ind(size(f_dists_ref_ref),...
    edges(1,:),edges(2,:)));


flow_lm = sample_with_flow(ref.xy, edges, source_idx, sink_idx, ...
    geo_dists, feat_dists, num_landmarks);

uniform_lm = sample_uniformly(ref.xy,num_landmarks);

f3 = figure; hold on;
scatter(ref.x, ref.y,'k','.', 'HandleVisibility','Off');
scatter(ref.x(flow_lm), ref.y(flow_lm),'o');
scatter(ref.x(uniform_lm), ref.y(uniform_lm),[],1:num_landmarks,'*');
legend({'Flow landmarks','Uniform landmarks'})
title('Selected landmarks')
xlabel('UTM Easting [m]')
ylabel('UTM Northing [m]')
hold off;
saveas(f3, fullfile('results', 'landmarks.eps'),'epsc')



%% Match with flow

F = f_dists_query_ref;
D = pdist2(query.xy', ref.xy');

% Match with flow

flow_matches = match_with_flow(query.xy,ref.xy(:,flow_lm),F(:,flow_lm));

% Retrieve without matching
[~, feature_matches] = sort(F(:,uniform_lm),2);


%% Display accuracy
f4 = figure; hold on;
plotAccVsDist(D(:,uniform_lm), feature_matches(:, 1:10), [0, 0.4470, 0.7410], ':');
plotAccVsDist(D(:,flow_lm), flow_matches, [0.4660, 0.6740, 0.1880], '-');
plotAccVsDist(D(:,uniform_lm), feature_matches(:, 1:1), [0.8500, 0.3250, 0.0980], ':');
legend({'Top-10 initial','Ours','Top-1 initial'})
title('Accuracy vs. distance')
xlabel('Distance [m]')
ylabel('Accuracy [%]')
grid on;
hold off;
saveas(f4, fullfile('results', 'accuracy.eps'),'epsc')



