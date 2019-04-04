function [psum] = integrate_path(xy)
%INTEGRATE_PATH(xy) where xy is a 2xn matrix of x and y coordinates

dx = sqrt(sum((xy(:, 2:end) - xy(:, 1:(end -1))).^2, 1));
psum = [0, dx];

for i = 2:length(psum)
    psum(i) = psum(i) + psum(i-1);
end

end
