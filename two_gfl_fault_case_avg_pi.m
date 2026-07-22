function two_gfl_fault_case_avg_pi()
%TWO_GFL_FAULT_CASE_AVG_PI
% Two grid-following inverters with a temporary short-circuit at bus 2.
% This version uses a reduced-order averaged current-control model:
%   PLL + PI current controller + current plant + network algebraic solve.
%
% Compared with the earlier version, the controller uses only partial
% terminal-voltage feedforward, so a fault can produce a visible current
% transient instead of being algebraically canceled.

clearvars -except nargin; %#ok<CLVAR>
close all;
clc;

p = paper_two_gfl_common_avg_pi();
caseCfg.type1 = 'GFL';
caseCfg.type2 = 'GFL';

% Initial guess.
x0 = struct();
x0.th1  = 0;
x0.xi1  = 0;
x0.Id1  = p.Idref1;
x0.Iq1  = p.Iqref1;
x0.zId1 = 0;
x0.zIq1 = 0;
x0.th2  = 0;
x0.xi2  = 0;
x0.Id2  = p.Idref2;
x0.Iq2  = p.Iqref2;
x0.zId2 = 0;
x0.zIq2 = 0;

% Align the initial state to the pre-fault equilibrium.
x0 = align_initial_state(x0, p, caseCfg);
xInit = [x0.th1; x0.xi1; x0.Id1; x0.Iq1; x0.zId1; x0.zIq1; ...
         x0.th2; x0.xi2; x0.Id2; x0.Iq2; x0.zId2; x0.zIq2];

opts = odeset('RelTol',1e-6,'AbsTol',1e-8,'MaxStep',1e-3);

[t1, x1] = ode15s(@(t,x) ode_case(t, x, p, caseCfg), [0 p.tFaultOn], xInit, opts);
xFaultOn = x1(end,:).';

[t2, x2] = ode15s(@(t,x) ode_case(t, x, p, caseCfg), [p.tFaultOn p.tFaultOff], xFaultOn, opts);
xFaultOff = x2(end,:).';

[t3, x3] = ode15s(@(t,x) ode_case(t, x, p, caseCfg), [p.tFaultOff p.tEnd], xFaultOff, opts);

% Stitch together, removing duplicate boundary points.
t = [t1; t2(2:end); t3(2:end)];
x = [x1; x2(2:end,:); x3(2:end,:)];

% Post-process.
N = numel(t);
store = init_store(N);
thetaDiff = zeros(N,1);

for k = 1:N
    s = state_to_struct(x(k,:).');
    faultActive = (t(k) >= p.tFaultOn) && (t(k) < p.tFaultOff);
    out = solve_two_bus_network_gfl_avg_pi(s, p, caseCfg, faultActive);

    store = push_store(store, k, out, s, p);
    thetaDiff(k) = wrap_pi(s.th1 - s.th2);
end

plot_case(t, store, thetaDiff, 'Two GFLs with temporary bus-2 fault (averaged PI current control, partial feedforward)');

end

function dx = ode_case(t, x, p, caseCfg)
    s = state_to_struct(x);
    faultActive = (t >= p.tFaultOn) && (t < p.tFaultOff);
    out = solve_two_bus_network_gfl_avg_pi(s, p, caseCfg, faultActive);

    % Inverter 1
    v1dq = out.V1 * exp(-1j * s.th1);
    v1d = real(v1dq);
    v1q = imag(v1dq);
    e1 = v1q;

    w1   = p.w0 + p.kppll * e1 + s.xi1;  % rad/s
    w1pu = w1 / p.w0;

    ed1 = p.Idref1 - s.Id1;
    eq1 = p.Iqref1 - s.Iq1;

    [vcd1, vcq1, awd1, awq1] = current_controller(v1d, v1q, s.Id1, s.Iq1, ed1, eq1, w1pu, p, s.zId1, s.zIq1);
    dzId1 = p.kic * ed1 + p.kaw * awd1;
    dzIq1 = p.kic * eq1 + p.kaw * awq1;

    dId1 = (vcd1 - v1d - p.Rf * s.Id1 + w1pu * p.Xf * s.Iq1) / p.Tlf;
    dIq1 = (vcq1 - v1q - p.Rf * s.Iq1 - w1pu * p.Xf * s.Id1) / p.Tlf;

    dxi1 = p.kipll * e1;
    dth1 = w1;

    % Inverter 2
    v2dq = out.V2 * exp(-1j * s.th2);
    v2d = real(v2dq);
    v2q = imag(v2dq);
    e2 = v2q;

    w2   = p.w0 + p.kppll * e2 + s.xi2;
    w2pu = w2 / p.w0;

    ed2 = p.Idref2 - s.Id2;
    eq2 = p.Iqref2 - s.Iq2;

    [vcd2, vcq2, awd2, awq2] = current_controller(v2d, v2q, s.Id2, s.Iq2, ed2, eq2, w2pu, p, s.zId2, s.zIq2);
    dzId2 = p.kic * ed2 + p.kaw * awd2;
    dzIq2 = p.kic * eq2 + p.kaw * awq2;

    dId2 = (vcd2 - v2d - p.Rf * s.Id2 + w2pu * p.Xf * s.Iq2) / p.Tlf;
    dIq2 = (vcq2 - v2q - p.Rf * s.Iq2 - w2pu * p.Xf * s.Id2) / p.Tlf;

    dxi2 = p.kipll * e2;
    dth2 = w2;

    dx = [dth1; dxi1; dId1; dIq1; dzId1; dzIq1; dth2; dxi2; dId2; dIq2; dzId2; dzIq2];
end

function [vcd, vcq, awd, awq] = current_controller(vd, vq, id, iq, ed, eq, wpu, p, zId, zIq)
    % PI current controller with decoupling/feedforward and vector saturation.
    % A partial voltage feedforward keeps the controller physical without
    % canceling the fault entirely.
    vcd_unsat = p.kff * vd + p.Rf * id - wpu * p.Xf * iq + p.kpc * ed + zId;
    vcq_unsat = p.kff * vq + p.Rf * iq + wpu * p.Xf * id + p.kpc * eq + zIq;

    [vcd, vcq] = saturate_vector(vcd_unsat, vcq_unsat, p.vConvMax);

    awd = vcd - vcd_unsat;
    awq = vcq - vcq_unsat;
end

function [ud, uq] = saturate_vector(ud, uq, umax)
    mag = hypot(ud, uq);
    if mag > umax
        s = umax / mag;
        ud = s * ud;
        uq = s * uq;
    end
end

function x0 = align_initial_state(x0, p, caseCfg)
    % Fixed-point alignment to the pre-fault equilibrium.
    for k = 1:40 %#ok<NASGU>
        s = state_to_struct([x0.th1; x0.xi1; x0.Id1; x0.Iq1; x0.zId1; x0.zIq1; ...
                             x0.th2; x0.xi2; x0.Id2; x0.Iq2; x0.zId2; x0.zIq2]);
        out = solve_two_bus_network_gfl_avg_pi(s, p, caseCfg, false);

        % Align PLL frames to the local bus voltages.
        x0.th1 = angle(out.V1);
        x0.th2 = angle(out.V2);

        % Set the current states to the commanded references.
        x0.Id1 = p.Idref1;
        x0.Iq1 = p.Iqref1;
        x0.Id2 = p.Idref2;
        x0.Iq2 = p.Iqref2;

        % Set current-controller integrators to the steady-state bias.
        v1dq = out.V1 * exp(-1j * x0.th1);
        v2dq = out.V2 * exp(-1j * x0.th2);
        v1d = real(v1dq);
        v1q = imag(v1dq);
        v2d = real(v2dq);
        v2q = imag(v2dq);

        x0.zId1 = (1 - p.kff) * v1d;
        x0.zIq1 = (1 - p.kff) * v1q;
        x0.zId2 = (1 - p.kff) * v2d;
        x0.zIq2 = (1 - p.kff) * v2q;

        % Set PLL integrators so omega starts nominal.
        e1 = v1q;
        e2 = v2q;
        x0.xi1 = -p.kppll * e1;
        x0.xi2 = -p.kppll * e2;
    end
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
    s.xi2  = x(8);
    s.Id2  = x(9);
    s.Iq2  = x(10);
    s.zId2 = x(11);
    s.zIq2 = x(12);

    % Kept for compatibility with older helper structures.
    s.E1d = NaN; s.E1q = NaN;
    s.E2d = NaN; s.E2q = NaN;
end

function store = init_store(N)
    store.w1  = zeros(N,1);
    store.w2  = zeros(N,1);
    store.id1 = zeros(N,1); store.iq1 = zeros(N,1);
    store.id2 = zeros(N,1); store.iq2 = zeros(N,1);
    store.vd1 = zeros(N,1); store.vq1 = zeros(N,1);
    store.vd2 = zeros(N,1); store.vq2 = zeros(N,1);
    store.P1  = zeros(N,1); store.P2  = zeros(N,1);
end

function store = push_store(store, k, out, s, p)
    e1 = imag(out.V1 * exp(-1j * s.th1));
    e2 = imag(out.V2 * exp(-1j * s.th2));

    w1 = p.w0 + p.kppll * e1 + s.xi1;
    w2 = p.w0 + p.kppll * e2 + s.xi2;

    store.w1(k) = w1 / p.w0;
    store.w2(k) = w2 / p.w0;

    store.id1(k) = out.id1; store.iq1(k) = out.iq1;
    store.id2(k) = out.id2; store.iq2(k) = out.iq2;
    store.vd1(k) = out.vd1; store.vq1(k) = out.vq1;
    store.vd2(k) = out.vd2; store.vq2(k) = out.vq2;
    store.P1(k)  = out.P1;  store.P2(k)  = out.P2;
end

function plot_case(t, store, thetaDiff, ttl)
    figure('Color','w','Name',ttl);
    tiledlayout(4,2,'TileSpacing','compact','Padding','compact');
    sgtitle(ttl);

    nexttile([1 2]);
    plot(t, store.w1, 'LineWidth', 1.1); hold on;
    plot(t, store.w2, 'LineWidth', 1.1);
    grid on;
    xline(0.20,'k--'); xline(0.26,'k--');
    title('PLL frequency');
    ylabel('\omega (pu)');
    legend('\omega_1','\omega_2','Location','best');

    nexttile;
    plot(t, store.id1, 'LineWidth', 1.1); hold on;
    plot(t, store.iq1, 'LineWidth', 1.1);
    grid on;
    xline(0.20,'k--'); xline(0.26,'k--');
    title('GFL1 current');
    legend('i_{d1}','i_{q1}','Location','best');

    nexttile;
    plot(t, store.id2, 'LineWidth', 1.1); hold on;
    plot(t, store.iq2, 'LineWidth', 1.1);
    grid on;
    xline(0.20,'k--'); xline(0.26,'k--');
    title('GFL2 current');
    legend('i_{d2}','i_{q2}','Location','best');

    nexttile;
    plot(t, store.vd1, 'LineWidth', 1.1); hold on;
    plot(t, store.vq1, 'LineWidth', 1.1);
    grid on;
    xline(0.20,'k--'); xline(0.26,'k--');
    title('GFL1 voltage');
    legend('v_{d1}','v_{q1}','Location','best');

    nexttile;
    plot(t, store.vd2, 'LineWidth', 1.1); hold on;
    plot(t, store.vq2, 'LineWidth', 1.1);
    grid on;
    xline(0.20,'k--'); xline(0.26,'k--');
    title('GFL2 voltage');
    legend('v_{d2}','v_{q2}','Location','best');

    nexttile([1 2]);
    plot(t, rad2deg(thetaDiff), 'LineWidth', 1.2);
    grid on;
    xline(0.20,'k--'); xline(0.26,'k--');
    xlabel('Time (s)');
    ylabel('\theta_1 - \theta_2 (deg)');
    title('Current-angle difference');
end

function a = wrap_pi(x)
    a = atan2(sin(x), cos(x));
end
