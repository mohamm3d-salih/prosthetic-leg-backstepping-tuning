function [xd, dxd, ddxd] = b2_ref_eval(t, R)
%B2_REF_EVAL  Shank reference and its first two derivatives at time t.
if R.type == 0
    xd = R.off + R.A*sin(R.w*t);
    dxd = R.A*R.w*cos(R.w*t);
    ddxd = -R.A*R.w^2*sin(R.w*t);
else
    W = 2*pi/R.period;
    k = 1:numel(R.a);
    c = cos(k*W*t);  s = sin(k*W*t);
    xd = R.a0 + sum(R.a.*c + R.b.*s);
    dxd = sum(k*W.*(-R.a.*s + R.b.*c));
    ddxd = -sum((k*W).^2.*(R.a.*c + R.b.*s));
end
end
