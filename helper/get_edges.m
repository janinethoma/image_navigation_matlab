function [edges, dists] = get_edges(xy, max_dist, diff_to_path)
%GET_EDGES Get edges with maximal length max_dist from 2xn node coordinate array
% This method is used for sequences where a floor plan is not available.

n = size(xy,2);

%create KD tree for the data xy'
tree = KDTreeSearcher(xy');

%find all the points in tree 
%that are within distance max_dist of the points in xy'
[Idx,D] = rangesearch(tree,xy',max_dist); % smaller or equals max_dist

% Count edges
%in Matlab 2016b, the length of D{i} is all equal to size(tree.X, 1) 
%that is, it can be viewd as a complete graph
%where the number of edges is (#node)*(#node-1)/2
count = 0;
for i = 1:n
    count = count + length(D{i}) -1;
end
num_edges = count/2;

%initialize edge and distance information
edges = NaN(2,num_edges);
dists = NaN(1,num_edges);

% Generate edges
%the variable "e" is the index for edges and dists
e = 1;
for i = 1:n
    for j = 2:length(Idx{i})
        %no edge generated between point i and point j
        if D{i}(j) == 0
            continue
        end
        
        %edge generated case
        %here, the conditional sentecne cannot be set as "else ..."
        %because the edge(A,B) is exactly the same as edge(B,A)
        if i < Idx{i}(j)
            edges(:,e) = [i;Idx{i}(j)];
            dists(e) = D{i}(j);
            e = e + 1;
        end
    end
end
edges(:,e:end) = [];
dists(e:end) = [];

%prepare the integrated path "l"
l = integrate_path(xy);

%if the (integrated path between the point A and point B - dist(A,B)) > diff_to_path
%cut this edge, since it will cut the corner
%for example,
%     ---------
%     |       |
%     A       |
%     |       B
% -----       -------
%as illustrated, edge(A,B) will cut the corner
to_delete = ((l(edges(2,:))-l(edges(1,:)))-dists)>diff_to_path;
disp(['Deleting ' num2str(sum(to_delete)) ' edges, as they cut corners.'])
edges(:,to_delete) = [];
dists(to_delete) = [];

end

