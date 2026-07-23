function two_gfm_fault_case_p_droop()
% TWO_GFM_FAULT_CASE_P_DROOP
% Two grid-forming inverters with P-w droop control, connected by an
% inductive tie line and subjected to a temporary short-circuit at bus 2.
%
% The model is intentionally reduced-order:
%   - voltage-forming sources behind virtual inductance
%   - low-pass filtered active power
%   - P-w droop control
%   - algebraic network solution at each time step

clearvars -except nargin; %#ok<CLVAR>
close all;
clc;

p = paper_two_gfm_common();
caseCfg.type1 = 'GFM';
caseCfg.type2 = 'GFM';

% Initial guess for the pre-fault equilibrium.
x0 = struct();
x0.th1 = 0;
x0.p1f = p.Pref1;
x0.th2 = 0;
x0.p2f = p.Pref2;

x0 = align_initial_state(x0, p, caseCfg);
xInit = [x0.th1; x0.p1f; x0.th2; x0.p2f];

opts = odeset('RelTol', p.relTol, 'AbsTol', p.absTol, 'MaxStep', p.maxStep);

[t1, x1] = ode15s(@(t, x) ode_case(t, x, p, caseCfg), [0 p.tFaultOn], xInit, opts);
xFaultOn = x1(end,:).';

[t2, x2] = ode15s(@(t, x) ode_case(t, x, p, caseCfg), [p.tFaultOn p.tFaultOff], xFaultOn, opts);
xFaultOff = x2(end,:).';

[t3, x3] = ode15s(@(t, x) ode_case(t, x, p, caseCfg), [p.tFaultOff p.tEnd], xFaultOff, opts);

% Stitch together, removing duplicate boundary points.
t = [t1; t2(2:end); t3(2:end)];
x = [x1; x2(2:end,:); x3(2:end,:)];

N = numel(t);
store = init_store(N);
thetaDiff = zeros(N,1);

for k = 1:N
    s = state_to_struct(x(k,:).');
    faultActive = (t(k) >= p.tFaultOn) && (t(k) < p.tFaultOff);
    out = solve_two_gfm_network(s, p, caseCfg, faultActive);

    store = push_store(store, k, out, s, p);
    thetaDiff(k) = wrap_pi(s.th1 - s.th2);
end

plot_case(t, store, thetaDiff, 'Two GFMs with temporary bus-2 fault');

end

function dx = ode_case(t, x, p, caseCfg)
    s = state_to_struct(x);
    faultActive = (t >= p.tFaultOn) && (t < p.tFaultOff);
    out = solve_two_gfm_network(s, p, caseCfg, faultActive);

    % Power filters.
    dp1f = p.wf * (out.P1 - s.p1f);
    dp2f = p.wf * (out.P2 - s.p2f);

    % P-w droop in per-unit frequency form.
    w1_pu = 1 + p.mDroop1 * (p.Pref1 - s.p1f);
    w2_pu = 1 + p.mDroop2 * (p.Pref2 - s.p2f);

    dth1 = p.w0 * w1_pu;
    dth2 = p.w0 * w2_pu;

    dx = [dth1; dp1f; dth2; dp2f];
end

function x0 = align_initial_state(x0, p, caseCfg)
    % Find a pre-fault angle difference that matches the desired power flow.
    delta0 = find_prefault_delta(p, caseCfg);
    x0.th1 = 0.5 * delta0;
    x0.th2 = -0.5 * delta0;
    x0.p1f = p.Pref1;
    x0.p2f = p.Pref2;

    % One consistency check pass.
    s = state_to_struct([x0.th1; x0.p1f; x0.th2; x0.p2f]);
    out = solve_two_gfm_network(s, p, caseCfg, false); %#ok<NASGU>
end

function delta0 = find_prefault_delta(p, caseCfg)
    grid = linspace(-1.4, 1.4, 281);
    vals = zeros(size(grid));
    for k = 1:numel(grid)
        vals(k) = power_mismatch(grid(k), p, caseCfg);
    end

    idx = find(vals(1:end-1) .* vals(2:end) <= 0, 1, 'first');
    if isempty(idx)
        [~, j] = min(abs(vals));
        delta0 = grid(j);
        return;
    end

    delta0 = fzero(@(d) power_mismatch(d, p, caseCfg), [grid(idx), grid(idx+1)]);
end

function err = power_mismatch(delta, p, caseCfg)
    s = struct();
    s.th1 = 0.5 * delta;
    s.th2 = -0.5 * delta;
    out = solve_two_gfm_network(s, p, caseCfg, false);
    err = out.P1 - p.Pref1;
end

function s = state_to_struct(x)
    x = x(:);
    s.th1 = x(1);
    s.p1f = x(2);
    s.th2 = x(3);
    s.p2f = x(4);
end

function store = init_store(N)
    store.w1 = zeros(N,1);
    store.w2 = zeros(N,1);
    store.P1 = zeros(N,1);
    store.P2 = zeros(N,1);
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
    w1_pu = 1 + p.mDroop1 * (p.Pref1 - s.p1f);
    w2_pu = 1 + p.mDroop2 * (p.Pref2 - s.p2f);

    store.w1(k) = w1_pu;
    store.w2(k) = w2_pu;
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
    plot(t, store.w1, 'LineWidth', 1.15); hold on;
    plot(t, store.w2, 'LineWidth', 1.15);
    grid on;
    xline(0.20,'k--'); xline(0.26,'k--');
    title('P-w droop frequency');
    ylabel('\omega (pu)');
    legend('\omega_1','\omega_2','Location','best');

    nexttile([1 2]);
    plot(t, store.P1, 'LineWidth', 1.15); hold on;
    plot(t, store.P2, 'LineWidth', 1.15);
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
    title('GFM1 current');
    legend('i_{d1}','i_{q1}','Location','best');

    nexttile;
    plot(t, store.id2, 'LineWidth', 1.1); hold on;
    plot(t, store.iq2, 'LineWidth', 1.1);
    grid on;
    xline(0.20,'k--'); xline(0.26,'k--');
    title('GFM2 current');
    legend('i_{d2}','i_{q2}','Location','best');

    nexttile;
    plot(t, store.vd1, 'LineWidth', 1.1); hold on;
    plot(t, store.vq1, 'LineWidth', 1.1);
    grid on;
    xline(0.20,'k--'); xline(0.26,'k--');
    title('GFM1 voltage');
    legend('v_{d1}','v_{q1}','Location','best');

    nexttile;
    plot(t, store.vd2, 'LineWidth', 1.1); hold on;
    plot(t, store.vq2, 'LineWidth', 1.1);
    grid on;
    xline(0.20,'k--'); xline(0.26,'k--');
    title('GFM2 voltage');
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
