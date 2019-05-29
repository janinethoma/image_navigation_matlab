function [sens] = sensitivity(edges, f_dists)
% SENSITIVITY(edges, f_dists) where edges is a 2xn array containing
% DIRECTED edges

nodes = 1:max(edges(1,:));
% It is unlikely, that a node has no edges.
% For data with many geographic outliers use
% nodes = unique(edges(1, :));
% and adjust code below
denominator = NaN(1, length(nodes));
for n = nodes
    denominator(n) = sum(f_dists(edges(1, :) == n));
end
sens = 1 - f_dists ./ denominator(edges(1, :));
end
