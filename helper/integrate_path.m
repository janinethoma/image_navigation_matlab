function [psum] = integrate_path(xy)
%INTEGRATE_PATH(xy) where xy is a 2xn matrix of x and y coordinates
%the concept of path integration: https://en.wikipedia.org/wiki/Path_integration

%first, calculate distance of each two adjacent points
%for example, there are three points named X, Y, Z
%distances we need to calculate are dist(X,Y) and dist(Y,Z)
dx = sqrt(sum((xy(:, 2:end) - xy(:, 1:(end -1))).^2, 1));

%next, put all the distances we calculated into psum
%that is, dist(X,Y) and dist(Y,Z) 
%and the order of the elements in psum will be [0, dist(X,Y), dist(Y,Z)]
psum = [0, dx];

%path integration
%the current term in psum will be the sum of current term and previous term
%that is, [0, (dist(X,Y) + 0), dist(Y,Z) + (dist(X,Y) + 0)]
for i = 2:length(psum)
    psum(i) = psum(i) + psum(i-1);
end

end
