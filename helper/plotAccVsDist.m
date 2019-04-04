function [] = plotAccVsDist(dist_geom,idx_match,clr,sty)
if nargin<3
    clr = [0 0 1];
    sty = '-';
end
allDist = zeros(size(dist_geom,1),1);

for i = 1:size(dist_geom,1)
    idx_i =  idx_match(i,:);
    allDist(i) = min(dist_geom(i,idx_i));
end

distSorted = sort(allDist);

maxVal = max(distSorted);
intrVal = (maxVal/100):(maxVal/100):maxVal;
accVal = zeros(size(intrVal));
N = length(intrVal);

for i = 1:N
    accVal(i) = 100*sum(distSorted<=intrVal(i))/length(allDist);
end
plot(intrVal,accVal,'LineWidth',2,'color',clr,'LineStyle',sty); hold on;
end