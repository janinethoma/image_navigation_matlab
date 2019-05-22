function [idx] = sample_with_flow(xy, edges, source_idx, sink_idx, ...
    geo_dists, feat_dists, num_to_choose)
%%SAMPLE_WITH_FLOW(xy, edges, source_idx, sink_idx,  geo_dists, feat_dists, num_to_choose)
% Find num landmarks using network flow
%
% INPUT:
% xy: 2xn node locations
% edges: 2xm edges represented as pair of node indices
% source_idx: indices of source nodes
% sink_idx: indices of sink nodes
% geo_dists: 1xm geometric length of edges
% feat_dists: 1xm feat distance between vertices connected by edge
% num_to_choose: number of landmarks to select
%
% OUTPUT:
% idx: node indices of selected landmarks


T = 0.1; % Total flow
a_dist = 4; % Distance between anchors
a_nn = 5; % Number of images in anchor nbh
tg = 0.1; % Flow through anchor nbh

%% Get costs
feat_dists = feat_dists / median(feat_dists);
costs = 1 ./ (1e-6 + feat_dists);

%% Get anchors
[anchors, a_nbh] = greedy_anchors(xy, a_dist, a_nn);

%% Get capacities
caps = geo_dists;

% Source capacity
for i = 1:length(source_idx)
    caps(edges(1, :) == source_idx(i)) = T;
    caps(edges(2, :) == source_idx(i)) = T;
end

% Sink capacity
for i = 1:length(sink_idx)
    caps(edges(1, :) == sink_idx(i)) = T;
    caps(edges(2, :) == sink_idx(i)) = T;
end

%% Convert to directed graph
arc_cap = [caps, caps];
arc_base = [costs, costs];
arc_i = [edges(1, :), edges(2, :)];
arc_j = [edges(2, :), edges(1, :)];

n = size(xy, 2);
narcs = size(arc_i, 2);
disp(['The number of nodes is ', num2str(n)]);
disp(['The number of edges is ', num2str(narcs)]);

%% Get sensitivity
% This should be done with directed graph
arc_sens = sensitivity([arc_i; arc_j], [feat_dists, feat_dists]);

%% Mosek code
disp('Starting optimization')
import mosek.fusion.*;
M = Model('LandmarkSelectionModel');
f = M.variable('f', narcs, Domain.inRange(0, arc_cap)); % Flow per edge
z = M.variable('z', narcs, Domain.greaterThan(0)); % Additional cost term

% Set the objective:
M.objective('Minimize total cost', ObjectiveSense.Minimize, ...
    Expr.add(Expr.dot(arc_base, f), Expr.sum(z)));

% Flow conservation constraints
for idx = 1:n % For each node
    f_tot = 0; % Total flow
    select_i = find(arc_i == idx);
    if ~isempty(select_i) % Edges outgoing from node idx
        f_tot = Expr.sum(f.pick(select_i));
    end
    select_j = find(arc_j == idx);
    if ~isempty(select_j) % Edges incoming to node idx
        f_tot = Expr.sub(f_tot, Expr.sum(f.pick(select_j)));
    end
    
    if isempty(union(select_i, select_j))
        continue;
    end
    
    if ~isempty(intersect(idx, source_idx))
        M.constraint(f_tot, Domain.equalsTo(T));
    elseif ~isempty(intersect(idx, sink_idx))
        M.constraint(f_tot, Domain.equalsTo(-T*(length(source_idx) / length(sink_idx))));
    else
        M.constraint(f_tot, Domain.equalsTo(0));
    end
end

% Anchor constraint, i.e. geometric representation
for a = 1:length(anchors)
    all_adj_edges = [];
    for nn = a_nbh(a, :) % Each node in given anchor nbh
        all_adj_edges = [all_adj_edges, find(arc_i == nn), find(arc_j == nn)];
    end
    all_adj_edges = unique(all_adj_edges);
    M.constraint(Expr.sum(f.pick(all_adj_edges)), Domain.greaterThan(tg));
end

% Visual representation
% Rotated quadratic cone = 2*lhs1*lhs2>rhs^2
lhs1 = Expr.mul(0.5, Expr.sub(arc_cap, f));
lhs2 = Expr.mulElm(z, (1 ./ (arc_sens .* arc_cap)));
stack = Expr.hstack(lhs1, lhs2, f);
M.constraint(stack, Domain.inRotatedQCone().axis(2)); % Each row is in a rotated quadratic cone

disp('Set constraints. Solving.')

M.solve();
M.dispose();

flow = f.level();


%% Choose nodes with highest flow
for i = 1:length(xy)
    node_flow(i) = sum(flow(arc_i == i));
end

[~, idx] = sort(node_flow, 'descend');

idx = idx(1:num_to_choose); % Implicitly choose tau to get desired number of landmarks
idx = sort(idx); % Return sorted landmarks
end
