function a = huberLoss(v,t)
if abs(v)<=t
    a = (v^2)/2;
else 
    a = t*(abs(v)-t/2);
end
end