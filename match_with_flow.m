function [matchFinalIdx] = match_with_flow(query_xy, ref_xy, visualDist)
%MATCH_WITH_FLOW(query_xy, ref_xy,  visual_Dist)
% 2xn query_xy 2xm ref_xy nxm visual_Dist

topN = 10;
distBoundInit = 25;

t = mean(ref_xy, 2);
s = mean(sqrt(sum((ref_xy - t).^2, 1)));
distBoundNorm = distBoundInit / s;

normMat = sqrt(2) * [1 / s, 0; 0, 1 / s];
query_xy = normMat * (query_xy - t);
ref_xy = normMat * (ref_xy - t);

visualDist = visualDist / median(min(visualDist));

[distValu1, idx_init1] = sort(visualDist, 2); % visual distance
[distValu2, idx_init2] = sort(visualDist, 1); % visual distance

%% Initialization
Vxy = ref_xy;
numP = size(query_xy, 2);
numG = size(Vxy, 2);
T = numP;
matchLenght = min([numP, numG, topN]);
costThresh = median(distValu1(:, 1));

nodeIdx = 0;
arc_cap = [];
arc_i = [];
arc_j = [];
arc_base = [];

%% S to G, s=1
for i = 1:numG
    arc_i = [arc_i, 1];
    arc_j = [arc_j, i + 1];
    arc_base = [arc_base, 0]; % no cost
    arc_cap = [arc_cap, numP];
end
nodeIdx = nodeIdx + 1;

%% form G to P
oneDirStart = length(arc_i) + 1;
% % % connect oneside
for i = 1:numG
    for j = 1:min(matchLenght, numP)
        arc_i = [arc_i, nodeIdx + i];
        arc_j = [arc_j, nodeIdx + numG + idx_init2(j, i)];
        arc_base = [arc_base, huberLoss(distValu2(j, i), costThresh)]; %huber cost
        arc_cap = [arc_cap, numP];
    end
end
% connect otherside
for i = 1:numP
    for j = 1:min(matchLenght, numG)
        arc_i = [arc_i, nodeIdx + idx_init1(i, j)];
        arc_j = [arc_j, nodeIdx + numG + i];
        arc_base = [arc_base, huberLoss(distValu1(i, j), costThresh)]; %huber cost
        arc_cap = [arc_cap, numP];
    end
end
oneDirEnd = length(arc_i);
nodeIdx = nodeIdx + numG;

%% form P to T
for i = 1:numP
    arc_i = [arc_i, nodeIdx + i];
    arc_j = [arc_j, nodeIdx + numP + 1];
    arc_base = [arc_base, 0]; % no cost
    arc_cap = [arc_cap, 1]; % limited capacity
end

%% Optimization part
mV = idx_init1(:, 1:matchLenght);
numAllV = 1 + numG + numP + 1; % 1 source, 1 sink, numG, numP
% Mosek code
narcs = size(arc_j, 2);
import mosek.fusion.*;
M = Model('You are amazing, just the way you are!');
x = M.variable('x', narcs, Domain.inRange(0, arc_cap));
y = M.variable('y', narcs, Domain.inRange(0, arc_cap));

% Set the objective:
M.objective('maximize your awesomeness', ...
    ObjectiveSense.Minimize, Expr.add(Expr.dot(arc_base, x), Expr.dot(arc_base, y)));

for idx = 1:numAllV
    v = 0;
    select_i = find(arc_i == idx);
    if ~isempty(select_i)
        v = Expr.sub(Expr.sum(x.pick(select_i)), Expr.sum(y.pick(select_i)));
    end
    select_j = find(arc_j == idx);
    if ~isempty(select_j)
        v = Expr.add(v, Expr.sub(Expr.sum(y.pick(select_j)), Expr.sum(x.pick(select_j))));
    end
    
    if isempty(union(select_i, select_j))
        continue;
    end
    
    if idx == 1 % source
        M.constraint(v, Domain.equalsTo(T));
    elseif idx == numAllV
        M.constraint(v, Domain.equalsTo(-T));
    else
        M.constraint(v, Domain.equalsTo(0));
    end
end

% Geometric constrains below
for i = numG + 2:numG + numP + 1
    selected1 = find(arc_j(oneDirStart:oneDirEnd) == i);
    selected2 = find(arc_j(oneDirStart:oneDirEnd) == i+1);
    
    if isempty(selected1) || isempty(selected2)
        continue;
    end
    flowIdx1 = selected1 + numG;
    flowIdx2 = selected2 + numG;
    
    vertexIdx1 = arc_i(selected1+numG) - 1;
    vertexIdx2 = arc_i(selected2+numG) - 1;
    
    v1All = Matrix.dense(Vxy(:, vertexIdx1));
    v2All = Matrix.dense(Vxy(:, vertexIdx2));
    
    x1All = x.pick(flowIdx1).transpose();
    x2All = x.pick(flowIdx2).transpose();
    
    fx1 = Expr.sum(Expr.mulElm(Expr.vstack(x1All, x1All), v1All), 1);
    fx2 = Expr.sum(Expr.mulElm(Expr.vstack(x2All, x2All), v2All), 1);
    
    M.constraint(Expr.vstack(distBoundNorm, Expr.sub(fx1, fx2)), Domain.inQCone());
end

%% Lets add additional constraints
jmpFlag = 1;
jumValue = 20;
if jmpFlag
    for i = numG + 2:numG + numP + 2 - jumValue
        selected1 = find(arc_j(oneDirStart:oneDirEnd) == i);
        selected2 = find(arc_j(oneDirStart:oneDirEnd) == i+jumValue);
        
        if isempty(selected1) || isempty(selected2)
            continue;
        end
        flowIdx1 = selected1 + numG;
        flowIdx2 = selected2 + numG;
        
        vertexIdx1 = arc_i(selected1+numG) - 1;
        vertexIdx2 = arc_i(selected2+numG) - 1;
        
        v1All = Matrix.dense(Vxy(:, vertexIdx1));
        v2All = Matrix.dense(Vxy(:, vertexIdx2));
        
        x1All = x.pick(flowIdx1).transpose();
        x2All = x.pick(flowIdx2).transpose();
        
        fx1 = Expr.sum(Expr.mulElm(Expr.vstack(x1All, x1All), v1All), 1);
        fx2 = Expr.sum(Expr.mulElm(Expr.vstack(x2All, x2All), v2All), 1);
        
        M.constraint(Expr.vstack(0.25*jumValue*distBoundNorm, Expr.sub(fx1, fx2)), Domain.inQCone());
    end
end

%%
M.solve();
M.dispose();

%%
flow = x.level();

%%  Computer matches from flow
matchFinalIdx = zeros(numP, 1);
for i = numG + 2:numG + numP + 1
    selected1 = find(arc_j(oneDirStart:oneDirEnd) == i);
    
    flowIdx1 = selected1 + numG;
    vertexIdx1 = arc_i(selected1+numG) - 1;
    
    v1All = Vxy(:, vertexIdx1);
    x1All = flow(flowIdx1);
    
    Y = sum(repmat(x1All', 2, 1).*v1All, 2); %% XL for error computation
    [~, idMin] = min(sum((ref_xy - Y).^2, 1)); %% idMin matched landmark
    matchFinalIdx(i-numG-1) = idMin;
end

end

