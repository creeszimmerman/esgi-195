
function x0 = align_gfl_gfm_exact(x0, p, caseCfg)
%ALIGN_GFL_GFM_EXACT
% Find a common pre-fault equilibrium for the mixed GFL-GFM system.
%
% We solve for the two source angles so that:
%   (1) the GFL q-axis terminal voltage is zero (PLL lock), and
%   (2) the GFM active power matches Pref2.
%
% After that, we initialize the GFL current-controller and PLL states, and
% set the GFM filtered power state.

th = [x0.th1; x0.th2];

for iter = 1:40 %#ok<NASGU>
    [res, out] = mixed_residual(th, p, caseCfg);

    if norm(res, inf) < 1e-11
        break;
    end

    J = zeros(2,2);
    for j = 1:2
        h = 1e-6 * (1 + abs(th(j)));
        thp = th;
        thp(j) = thp(j) + h;
        resp = mixed_residual(thp, p, caseCfg);
        J(:,j) = (resp - res) / h;
    end

    if rcond(J) < 1e-12
        step = -pinv(J) * res;
    else
        step = -J \ res;
    end

    % Conservative damping to avoid jumping to a different branch.
    stepNorm = norm(step, inf);
    if stepNorm > 0.5
        step = 0.5 * step / stepNorm;
    end

    th = th + step;
end

x0.th1 = th(1);
x0.th2 = th(2);

% Final network solve at the pre-fault equilibrium.
s = state_to_struct([x0.th1; x0.xi1; x0.Id1; x0.Iq1; x0.zId1; x0.zIq1; x0.th2; x0.p2f]);
out = solve_two_bus_network_gfl_gfm_exact(s, p, caseCfg, false);

% GFL current states and PLL.
x0.Id1 = p.Idref1;
x0.Iq1 = p.Iqref1;
v1dq = out.V1 * exp(-1j * x0.th1);
v1d = real(v1dq);
v1q = imag(v1dq);
x0.xi1 = -p.kppll * v1q;
x0.zId1 = (1 - p.kff) * v1d;
x0.zIq1 = (1 - p.kff) * v1q;

% GFM filtered power state.
x0.p2f = out.P2;
end

function [res, out] = mixed_residual(th, p, caseCfg)
    s = struct();
    s.th1 = th(1);
    s.th2 = th(2);
    s.Id1 = p.Idref1;
    s.Iq1 = p.Iqref1;

    out = solve_two_bus_network_gfl_gfm_exact(s, p, caseCfg, false);

    v1dq = out.V1 * exp(-1j * th(1));
    v1q = imag(v1dq);

    res = [v1q; out.P2 - p.Pref2];
end

function s = state_to_struct(x)
    x = x(:);
    s.th1  = x(1);
    s.xi1  = x(2);
    s.Id1  = x(3);
    s.Iq1  = x(4);
    s.zId1 = x(5);
    s.zIq1 = x(6);
    s.th2  = x(7);
    s.p2f  = x(8);
end
