function [edges, dists] = get_edges(xy, max_dist, diff_to_path)
%GET_EDGES Get edges with maximal length max_dist from 2xn node coordinate array
% This method is used for sequences where a floor plan is not available.

n = size(xy,2);
tree = KDTreeSearcher(xy');
[Idx,D] = rangesearch(tree,xy',max_dist); % smaller or equals max_dist

% Count edges
count = 0;
for i = 1:n
    count = count + length(D{i}) -1;
end
num_edges = count/2;

edges = NaN(2,num_edges);
dists = NaN(1,num_edges);

% Generate edges
e = 1;
for i = 1:n
    for j = 2:length(Idx{i})
        if D{i}(j) == 0
            continue
        end
        
        if i < Idx{i}(j)
            edges(:,e) = [i;Idx{i}(j)];
            dists(e) = D{i}(j);
            e = e + 1;
        end
    end
end
edges(:,e:end) = [];
dists(e:end) = [];

l = integrate_path(xy);
to_delete = ((l(edges(2,:))-l(edges(1,:)))-dists)>diff_to_path;
disp(['Deleting ' num2str(sum(to_delete)) ' edges, as they cut corners.'])
edges(:,to_delete) = [];
dists(to_delete) = [];

end

