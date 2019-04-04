function [sens] = sensitivity(edges, f_dists)
% SENSITIVITY(edges, f_dists) where edges is a 2xn array containing
% DIRECTED edges

nodes = unique(edges(1, :));
denominator = NaN(1, length(nodes));
for n = nodes
    denominator(n) = sum(f_dists(edges(1, :) == n));
end
sens = 1 - f_dists ./ denominator(edges(1, :));
end
