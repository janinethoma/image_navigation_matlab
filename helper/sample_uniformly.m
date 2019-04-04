function [idx] = sample_uniformly(xy, num_to_take)
% SAMPLE_UNIFORMLY(xy,num_to_take)
% xy: 2 x n

l = integrate_path([xy(1, :); xy(2, :)]);
locs = linspace(0, l(end), num_to_take);

idx = NaN(1, num_to_take);
for j = 1:num_to_take
    [~, idx(j)] = min(abs(l-locs(j)));
    
end

end
