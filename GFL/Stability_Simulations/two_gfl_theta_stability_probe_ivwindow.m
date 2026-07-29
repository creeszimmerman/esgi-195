function two_gfl_theta_stability_probe_ivwindow()
%TWO_GFL_THETA_STABILITY_PROBE_IVWINDOW
% Start the two-GFL averaged PI model from a chosen initial angle offset
% theta1-theta2 and observe whether the trajectory settles.
%
% No short-circuit is applied in this version.
% The script prints the initial theta and Q difference, then plots:
%   1) theta1-theta2 over time
%   2) Q1 and Q2 over time
%   3) Q1-Q2 over time
%   4) I-V locus over the last user-selected time window
%
% Edit deltaThetaDeg below to probe different initial conditions.
% Edit ivWindowSec below to choose how much of the tail to plot for the
% I-V locus.

clearvars -except nargin; %#ok<CLVAR>
close all;
clc;

p = paper_two_gfl_common_avg_pi_local();
caseCfg.type1 = 'GFL';
caseCfg.type2 = 'GFL';

% Initial angle perturbation to test.
deltaThetaDeg = 100;
deltaTheta = deg2rad(deltaThetaDeg);

% Time window for the I-V plot at the end of the simulation.
ivWindowSec = 0.1;  % change this to whatever you want

% Base initial condition.
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

% First align to the pre-fault equilibrium.
x0 = align_initial_state(x0, p, caseCfg);

% Then apply a symmetric angle offset while keeping the mean angle fixed.
x0.th1 = x0.th1 + 0.5 * deltaTheta;
x0.th2 = x0.th2 - 0.5 * deltaTheta;

xInit = [x0.th1; x0.xi1; x0.Id1; x0.Iq1; x0.zId1; x0.zIq1; ...
         x0.th2; x0.xi2; x0.Id2; x0.Iq2; x0.zId2; x0.zIq2];

% Print the initial conditions.
out0 = solve_two_bus_network_gfl_avg_pi_local(state_to_struct(xInit), p, caseCfg);
thetaDiff0 = wrap_pi(x0.th1 - x0.th2);
Qdiff0 = out0.Q1 - out0.Q2;
fprintf('Initial theta1-theta2 = %.6f deg\n', rad2deg(thetaDiff0));
fprintf('Initial Q1-Q2         = %.6f pu\n', Qdiff0);
fprintf('Initial Q1            = %.6f pu\n', out0.Q1);
fprintf('Initial Q2            = %.6f pu\n', out0.Q2);

opts = odeset('RelTol',1e-6,'AbsTol',1e-8,'MaxStep',1e-3);
[t, x] = ode15s(@(t,x) ode_case(t, x, p, caseCfg), [0 p.tEnd], xInit, opts);

% Post-process.
N = numel(t);
store = init_store(N);
thetaDiff = zeros(N,1);
Qdiff = zeros(N,1);

for k = 1:N
    s = state_to_struct(x(k,:).');
    out = solve_two_bus_network_gfl_avg_pi_local(s, p, caseCfg);
    store = push_store(store, k, out, s, p);
    thetaDiff(k) = wrap_pi(s.th1 - s.th2);
    Qdiff(k) = out.Q1 - out.Q2;
end

% Simple boundedness screen (heuristic, not a formal proof).
stateMax = max([max(abs(thetaDiff)), max(abs(Qdiff))]);
if isfinite(stateMax) && stateMax < 1e3
    fprintf('Trajectory remained bounded over the simulation window.\n');
else
    fprintf('Trajectory grew very large or became non-finite.\n');
end

plot_case(t, store, thetaDiff, Qdiff, sprintf('Two GFLs: initial theta1-theta2 = %.1f deg', deltaThetaDeg));
plot_iv_window(t, store, ivWindowSec);
end

function dx = ode_case(~, x, p, caseCfg)
    s = state_to_struct(x);
    out = solve_two_bus_network_gfl_avg_pi_local(s, p, caseCfg);

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
    for iter = 1:40 %#ok<NASGU>
        s = state_to_struct([x0.th1; x0.xi1; x0.Id1; x0.Iq1; x0.zId1; x0.zIq1; ...
                             x0.th2; x0.xi2; x0.Id2; x0.Iq2; x0.zId2; x0.zIq2]);
        out = solve_two_bus_network_gfl_avg_pi_local(s, p, caseCfg);

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

function out = solve_two_bus_network_gfl_avg_pi_local(x, p, caseCfg)
    %SOLVE_TWO_BUS_NETWORK_GFL_AVG_PI_LOCAL
    th1 = x.th1;
    th2 = x.th2;

    switch caseCfg.type1
        case 'GFL'
            I1dq = x.Id1 + 1j * x.Iq1;
            Iinj1 = I1dq * exp(1j * th1);
        otherwise
            error('Only GFL-GFL is supported.');
    end

    switch caseCfg.type2
        case 'GFL'
            I2dq = x.Id2 + 1j * x.Iq2;
            Iinj2 = I2dq * exp(1j * th2);
        otherwise
            error('Only GFL-GFL is supported.');
    end

    % Two-bus nodal network:
    %   bus 1 -- Z12 -- bus 2
    %   |                 |
    %  Gload1           Gload2
    %   |                 |
    % ground            ground
    Y12 = 1 / p.Z12;
    Y1  = p.Gload1;
    Y2  = p.Gload2;

    Ybus = [Y1 + Y12, -Y12;
            -Y12,     Y2 + Y12];

    V = Ybus \ [Iinj1; Iinj2];
    V1 = V(1);
    V2 = V(2);

    % dq variables in each inverter's own rotating frame.
    v1dq = V1 * exp(-1j * th1);
    v2dq = V2 * exp(-1j * th2);
    i1dq = Iinj1 * exp(-1j * th1);
    i2dq = Iinj2 * exp(-1j * th2);

    out.V1 = V1; out.V2 = V2;
    out.I1 = Iinj1; out.I2 = Iinj2;
    out.vd1 = real(v1dq); out.vq1 = imag(v1dq);
    out.vd2 = real(v2dq); out.vq2 = imag(v2dq);
    out.id1 = real(i1dq); out.iq1 = imag(i1dq);
    out.id2 = real(i2dq); out.iq2 = imag(i2dq);
    out.P1 = real(V1 * conj(Iinj1));
    out.P2 = real(V2 * conj(Iinj2));
    out.Q1 = imag(V1 * conj(Iinj1));
    out.Q2 = imag(V2 * conj(Iinj2));
end

function store = init_store(N)
    store.w1  = zeros(N,1);
    store.w2  = zeros(N,1);
    store.id1 = zeros(N,1); store.iq1 = zeros(N,1);
    store.id2 = zeros(N,1); store.iq2 = zeros(N,1);
    store.vd1 = zeros(N,1); store.vq1 = zeros(N,1);
    store.vd2 = zeros(N,1); store.vq2 = zeros(N,1);
    store.Q1  = zeros(N,1); store.Q2  = zeros(N,1);
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
    store.Q1(k)  = out.Q1;  store.Q2(k)  = out.Q2;
end

function plot_case(t, store, thetaDiff, Qdiff, ttl)
    figure('Color','w','Name',ttl);
    tiledlayout(3,1,'TileSpacing','compact','Padding','compact');
    sgtitle(ttl);

    nexttile;
    plot(t, rad2deg(thetaDiff), 'LineWidth', 1.2);
    grid on;
    xlabel('Time (s)');
    ylabel('\theta_1 - \theta_2 (deg)');
    title('Initial angle difference response');

    nexttile;
    plot(t, store.Q1, 'LineWidth', 1.1); hold on;
    plot(t, store.Q2, 'LineWidth', 1.1);
    grid on;
    xlabel('Time (s)');
    ylabel('Q (pu)');
    legend('Q_1','Q_2','Location','best');
    title('Reactive power outputs');

    nexttile;
    plot(t, Qdiff, 'LineWidth', 1.2);
    grid on;
    xlabel('Time (s)');
    ylabel('Q_1 - Q_2 (pu)');
    title('Reactive power difference');
end

function plot_iv_window(t, store, windowSec)
    % Plot |I| vs |V| for the last user-selected time window.
    if ~isscalar(windowSec) || ~isfinite(windowSec) || windowSec <= 0
        error('ivWindowSec must be a positive finite scalar.');
    end

    t0 = max(t(1), t(end) - windowSec);
    idx = t >= t0;
    if nnz(idx) < 2
        idx = true(size(t));
        warning('Selected I-V window had fewer than 2 points; plotting the full simulation instead.');
    end

    V1mag = hypot(store.vd1(idx), store.vq1(idx));
    I1mag = hypot(store.id1(idx), store.iq1(idx));
    V2mag = hypot(store.vd2(idx), store.vq2(idx));
    I2mag = hypot(store.id2(idx), store.iq2(idx));
    tt = t(idx);

    figure('Color','w','Name',sprintf('I-V locus (last %.4f s)', windowSec));
    tiledlayout(1,2,'TileSpacing','compact','Padding','compact');

    nexttile;
    scatter(V1mag, I1mag, 18, tt, 'filled'); hold on;
    plot(V1mag, I1mag, 'k-', 'LineWidth', 1.1);
    plot(V1mag(1), I1mag(1), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 8);
    plot(V1mag(end), I1mag(end), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
    grid on;
    xlabel('|V_1| (pu)');
    ylabel('|I_1| (pu)');
    title(sprintf('GFL1 I-V locus (last %.4f s)', windowSec));
    colorbar;
    colormap(parula);
    legend('Samples','Trajectory','Initial','Final','Location','best');

    nexttile;
    scatter(V2mag, I2mag, 18, tt, 'filled'); hold on;
    plot(V2mag, I2mag, 'k-', 'LineWidth', 1.1);
    plot(V2mag(1), I2mag(1), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 8);
    plot(V2mag(end), I2mag(end), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
    grid on;
    xlabel('|V_2| (pu)');
    ylabel('|I_2| (pu)');
    title(sprintf('GFL2 I-V locus (last %.4f s)', windowSec));
    colorbar;
    colormap(parula);
    legend('Samples','Trajectory','Initial','Final','Location','best');
end

function a = wrap_pi(x)
    a = atan2(sin(x), cos(x));
end

function p = paper_two_gfl_common_avg_pi_local()
    % Parameters for a reduced-order two-GFL case with averaged current
    % control (PLL + PI current controller + current plant).

    p.f0 = 50;
    p.w0 = 2*pi*p.f0;

    % PLL parameters (same nominal values as the paper's Table II style setup).
    p.wpll  = 2*pi*15;
    p.kppll = p.wpll;
    p.kipll = p.wpll^2/4;

    % Current-control / output-filter surrogate.
    % Tlf is an effective current-plant time constant in seconds.
    p.Tlf = 0.018;
    p.Rf  = 0.03;
    p.Xf  = 0.12;

    % PI current controller gains (conservative values).
    p.kpc = 5;
    p.kic = 1000;

    % Partial voltage feedforward coefficient.
    p.kff = 0.0000001;

    % Simple back-calculation anti-windup gain.
    p.kaw = 6.0;

    % Optional converter-voltage saturation.
    p.vConvMax = 3;

    % Grid-following current references.
    p.Idref1 = 0.5;
    p.Iqref1 = 0.5; % 0.0;
    p.Idref2 = 0.5;
    p.Iqref2 = -0.5; %0.0;

    % Network model.
    p.Vref = 1.0;
    p.Gload1 = 0.4;
    p.Gload2 = 0.6;
    p.Z12    = 0.5 + 1j*0.3;

    % Optional overall current ceiling used only as a safety check.
    p.iMax = 5.0;

    % Simulation window.
    p.tEnd = 0.60;
end
