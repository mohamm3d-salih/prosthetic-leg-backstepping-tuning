function m = b2_metrics(out, R, ss_window)
%B2_METRICS  Performance measures of one simulation (degrees, rad.s, N.m).
%   RMSE      RMS shank-angle error over the whole run (deg)
%   RMSE_ss   RMS error over the last ss_window seconds (deg)
%   ess       max |error| over the last ss_window seconds (deg)
%   IAE, ITAE integral measures of the error (rad.s, rad.s^2)
%   umax, urms peak and RMS motor torque (N.m)
%   hip_pp    peak-to-peak hip angle over the last ss_window seconds (deg)
%   hip_max, hipvel_max  peak |hip angle| (deg) and |hip velocity| (deg/s)
%   flex_min, flex_max   min/max knee flexion th1 - th2 (deg); negative = hyperextension
if nargin < 3, ss_window = 10; end
xd = arrayfun(@(tt) b2_ref_eval(tt, R), out.t);
e = out.x(:, 3) - xd;  ed = e*180/pi;
last = out.t >= out.t(end) - ss_window;
flex = (out.x(:, 1) - out.x(:, 3))*180/pi;
m = struct('RMSE', sqrt(mean(ed.^2)), 'RMSE_ss', sqrt(mean(ed(last).^2)), ...
           'ess', max(abs(ed(last))), 'IAE', trapz(out.t, abs(e)), ...
           'ITAE', trapz(out.t, out.t.*abs(e)), 'umax', max(abs(out.u)), ...
           'urms', sqrt(mean(out.u.^2)), ...
           'hip_pp', (max(out.x(last, 1)) - min(out.x(last, 1)))*180/pi, ...
           'hip_max', max(abs(out.x(:, 1)))*180/pi, ...
           'hipvel_max', max(abs(out.x(:, 2)))*180/pi, ...
           'flex_min', min(flex), 'flex_max', max(flex), 'ok', out.ok);
end
