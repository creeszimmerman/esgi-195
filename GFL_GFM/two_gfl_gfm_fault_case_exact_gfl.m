
function two_gfl_gfm_fault_case_exact_gfl()
%TWO_GFL_GFM_FAULT_CASE_EXACT_GFL
% Mixed GFL-GFM case with an exact reduced-order GFL model and an exact
% reduced-order GFM model. A temporary short-circuit is applied at bus 2.

clearvars -except nargin; %#ok<CLVAR>
close all;
clc;

p = paper_gfl_gfm_exact_common();
caseCfg.type1 = 'GFL';
caseCfg.type2 = 'GFM';

% Initial guess.
x0 = struct();
x0.th1  = 0;
x0.xi1  = 0;
x0.Id1  = p.Idref1;
x0.Iq1  = p.Iqref1;
x0.zId1 = 0;
x0.zIq1 = 0;
x0.th2  = 0;
x0.p2f  = p.Pref2;

% Align both devices to the same pre-fault equilibrium.
x0 = align_gfl_gfm_exact(x0, p, caseCfg);
xInit = [x0.th1; x0.xi1; x0.Id1; x0.Iq1; x0.zId1; x0.zIq1; x0.th2; x0.p2f];

opts = odeset('RelTol',p.relTol,'AbsTol',p.absTol,'MaxStep',p.maxStep);

[t1, x1] = ode15s(@(t,x) ode_case(t, x, p, caseCfg, false), [0 p.tFaultOn], xInit, opts);
xFaultOn = x1(end,:).';

[t2, x2] = ode15s(@(t,x) ode_case(t, x, p, caseCfg, true), [p.tFaultOn p.tFaultOff], xFaultOn, opts);
xFaultOff = x2(end,:).';

[t3, x3] = ode15s(@(t,x) ode_case(t, x, p, caseCfg, false), [p.tFaultOff p.tEndFault], xFaultOff, opts);

% Stitch together, removing duplicate boundary points.
t = [t1; t2(2:end); t3(2:end)];
x = [x1; x2(2:end,:); x3(2:end,:)];

N = numel(t);
store = init_store(N);
thetaDiff = zeros(N,1);

for k = 1:N
    s = state_to_struct(x(k,:).');
    faultActive = (t(k) >= p.tFaultOn) && (t(k) < p.tFaultOff);
    out = solve_two_bus_network_gfl_gfm_exact(s, p, caseCfg, faultActive);

    store = push_store(store, k, out, s, p);
    thetaDiff(k) = wrap_pi(s.th1 - s.th2);
end

plot_case(t, store, thetaDiff, 'Mixed GFL-GFM with temporary bus-2 fault (exact GFL physics)');
end

function dx = ode_case(t, x, p, caseCfg, faultActive)
    s = state_to_struct(x);
    out = solve_two_bus_network_gfl_gfm_exact(s, p, caseCfg, faultActive);

    % --- GFL at bus 1 ---
    v1dq = out.V1 * exp(-1j * s.th1);
    v1d = real(v1dq);
    v1q = imag(v1dq);
    e1 = v1q;

    w1   = p.w0 + p.kppll * e1 + s.xi1;
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

    % --- GFM at bus 2 ---
    dp2f = p.wf * (out.P2 - s.p2f);
    w2_pu = 1 + p.mDroop2 * (p.Pref2 - s.p2f);
    dth2 = p.w0 * w2_pu;

    dx = [dth1; dxi1; dId1; dIq1; dzId1; dzIq1; dth2; dp2f];
end

function [vcd, vcq, awd, awq] = current_controller(vd, vq, id, iq, ed, eq, wpu, p, zId, zIq)
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

function store = init_store(N)
    store.w1  = zeros(N,1);
    store.w2  = zeros(N,1);
    store.P1  = zeros(N,1);
    store.P2  = zeros(N,1);
    store.vd1 = zeros(N,1);
    store.vq1 = zeros(N,1);
    store.vd2 = zeros(N,1);
    store.vq2 = zeros(N,1);
    store.id1 = zeros(N,1);
    store.iq1 = zeros(N,1);
    store.id2 = zeros(N,1);
    store.iq2 = zeros(N,1);
    store.V1mag = zeros(N,1);
    store.V2mag = zeros(N,1);
    store.I1mag = zeros(N,1);
    store.I2mag = zeros(N,1);
end

function store = push_store(store, k, out, s, p)
    e1 = imag(out.V1 * exp(-1j * s.th1));
    w1 = p.w0 + p.kppll * e1 + s.xi1;
    w2 = p.w0 + p.mDroop2 * 0; %#ok<NASGU>

    store.w1(k) = w1 / p.w0;
    store.w2(k) = (p.w0 + p.mDroop2 * (p.Pref2 - s.p2f)) / p.w0;
    store.P1(k) = out.P1;
    store.P2(k) = out.P2;
    store.vd1(k) = out.vd1;
    store.vq1(k) = out.vq1;
    store.vd2(k) = out.vd2;
    store.vq2(k) = out.vq2;
    store.id1(k) = out.id1;
    store.iq1(k) = out.iq1;
    store.id2(k) = out.id2;
    store.iq2(k) = out.iq2;
    store.V1mag(k) = out.V1mag;
    store.V2mag(k) = out.V2mag;
    store.I1mag(k) = out.I1mag;
    store.I2mag(k) = out.I2mag;
end

function plot_case(t, store, thetaDiff, ttl)
    figure('Color','w','Name',ttl);
    tiledlayout(5,2,'TileSpacing','compact','Padding','compact');
    sgtitle(ttl);

    nexttile([1 2]);
    plot(t, store.w1, 'LineWidth', 1.1); hold on;
    plot(t, store.w2, 'LineWidth', 1.1);
    grid on;
    xline(0.20,'k--'); xline(0.26,'k--');
    title('PLL frequency / droop frequency');
    ylabel('\omega (pu)');
    legend('\omega_1','\omega_2','Location','best');

    nexttile([1 2]);
    plot(t, store.P1, 'LineWidth', 1.1); hold on;
    plot(t, store.P2, 'LineWidth', 1.1);
    grid on;
    xline(0.20,'k--'); xline(0.26,'k--');
    title('Active power');
    ylabel('P (pu)');
    legend('P_1','P_2','Location','best');

    nexttile;
    plot(t, store.id1, 'LineWidth', 1.1); hold on;
    plot(t, store.iq1, 'LineWidth', 1.1);
    grid on;
    xline(0.20,'k--'); xline(0.26,'k--');
    title('GFL current');
    legend('i_{d1}','i_{q1}','Location','best');

    nexttile;
    plot(t, store.id2, 'LineWidth', 1.1); hold on;
    plot(t, store.iq2, 'LineWidth', 1.1);
    grid on;
    xline(0.20,'k--'); xline(0.26,'k--');
    title('GFM current');
    legend('i_{d2}','i_{q2}','Location','best');

    nexttile;
    plot(t, store.vd1, 'LineWidth', 1.1); hold on;
    plot(t, store.vq1, 'LineWidth', 1.1);
    grid on;
    xline(0.20,'k--'); xline(0.26,'k--');
    title('GFL voltage');
    legend('v_{d1}','v_{q1}','Location','best');

    nexttile;
    plot(t, store.vd2, 'LineWidth', 1.1); hold on;
    plot(t, store.vq2, 'LineWidth', 1.1);
    grid on;
    xline(0.20,'k--'); xline(0.26,'k--');
    title('GFM voltage');
    legend('v_{d2}','v_{q2}','Location','best');

    nexttile([1 2]);
    plot(t, rad2deg(thetaDiff), 'LineWidth', 1.2);
    grid on;
    xline(0.20,'k--'); xline(0.26,'k--');
    xlabel('Time (s)');
    ylabel('\theta_1 - \theta_2 (deg)');
    title('Angle difference');
end

function a = wrap_pi(x)
    a = atan2(sin(x), cos(x));
end
