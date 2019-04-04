function [selected, a_nbh] = greedy_anchors(xy, a_dist, a_nn)
% GREEDY_ANCHORS(xy, a_dist, a_nn) where xy has row dimensions and column samples

selected = randi(length(xy), 1, 1);
remaining = 1:length(xy);
remaining(remaining == selected) = [];

while length(selected) < length(xy)
    
    tree = KDTreeSearcher(xy(:, selected)');
    [~, D] = knnsearch(tree, xy(:, remaining)');
    
    [m, i] = max(D);
    
    if m < a_dist / 2
        break
    end
    
    selected = [selected, remaining(i)];
    remaining(remaining == remaining(i)) = [];
    
end

selected = sort(selected);
tree = KDTreeSearcher(xy(:, remaining)');

a_nbh = NaN(length(selected), 1+a_nn);

for a = 1:length(selected)
    a_nbh(a, :) = [selected(a), remaining(knnsearch(tree, xy(:, selected(a))', 'K', a_nn))];
    
end

end
