function two_gfl_theta_stability_probe()
%TWO_GFL_THETA_STABILITY_PROBE
% Standalone no-fault stability probe for the reduced two-GFL model.
% The script:
%   1) builds the nominal parameter set,
%   2) finds the pre-perturbation equilibrium,
%   3) applies a chosen initial theta1-theta2 offset,
%   4) simulates the no-fault response,
%   5) prints the initial theta difference and Q1-Q2,
%   6) plots time traces and the phase portrait Delta(theta) vs Delta(Q).

clearvars -except nargin; %#ok<CLVAR>
close all;
clc;

if exist('paper_two_gfl_common_avg_pi', 'file') ~= 2
    error('Could not find paper_two_gfl_common_avg_pi.m on the MATLAB path.');
end

p = paper_two_gfl_common_avg_pi();
caseCfg.type1 = 'GFL';
caseCfg.type2 = 'GFL';

% For this stability probe, keep the finite loads used in the more physical case.
p.Gload1 = 0.4;
p.Gload2 = 0.6;

% Choose the initial perturbation in theta1-theta2 (degrees).
thetaOffsetDeg = 100;  % edit this value as needed
thetaOffsetRad = deg2rad(thetaOffsetDeg);

% Initial guess for the equilibrium search.
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

% Align to the no-perturbation equilibrium first.
x0 = align_initial_state(x0, p, caseCfg);

% Apply the desired initial angle difference while keeping the average angle fixed.
x0.th1 = x0.th1 + thetaOffsetRad/2;
x0.th2 = x0.th2 - thetaOffsetRad/2;

xInit = [x0.th1; x0.xi1; x0.Id1; x0.Iq1; x0.zId1; x0.zIq1; ...
         x0.th2; x0.xi2; x0.Id2; x0.Iq2; x0.zId2; x0.zIq2];

% Print the initial condition summary.
sInit = state_to_struct(xInit);
outInit = solve_two_bus_network_gfl_avg_pi(sInit, p, caseCfg, false);
thetaDiff0 = wrap_pi(sInit.th1 - sInit.th2);
Qdiff0 = outInit.Q1 - outInit.Q2;

fprintf('\n=== Initial condition summary ===\n');
fprintf('Requested theta1-theta2 offset : %+8.3f deg\n', thetaOffsetDeg);
fprintf('Actual theta1-theta2 at t = 0  : %+8.3f deg\n', rad2deg(thetaDiff0));
fprintf('Initial Q1                     : %+12.6f pu\n', outInit.Q1);
fprintf('Initial Q2                     : %+12.6f pu\n', outInit.Q2);
fprintf('Initial Q1-Q2                  : %+12.6f pu\n', Qdiff0);
fprintf('=================================\n\n');

opts = odeset('RelTol',1e-6,'AbsTol',1e-8,'MaxStep',1e-3);

[t, x] = ode15s(@(tt,xx) ode_case(tt, xx, p, caseCfg), [0 p.tEnd], xInit, opts);

% Post-process.
N = numel(t);
store = init_store(N);
thetaDiff = zeros(N,1);
Qdiff = zeros(N,1);
deltaOmega = zeros(N,1);

for k = 1:N
    s = state_to_struct(x(k,:).');
    out = solve_two_bus_network_gfl_avg_pi(s, p, caseCfg, false);

    store = push_store(store, k, out, s, p);

    thetaDiff(k) = wrap_pi(s.th1 - s.th2);
    Qdiff(k) = out.Q1 - out.Q2;
    deltaOmega(k) = store.w1(k) - store.w2(k);
end

thetaDiffDeg = rad2deg(thetaDiff);

% Final values.
thetaDiffF = thetaDiffDeg(end);
QdiffF = Qdiff(end);

fprintf('=== Final condition summary ===\n');
fprintf('Final theta1-theta2            : %+8.3f deg\n', thetaDiffF);
fprintf('Final Q1                        : %+12.6f pu\n', store.Q1(end));
fprintf('Final Q2                        : %+12.6f pu\n', store.Q2(end));
fprintf('Final Q1-Q2                     : %+12.6f pu\n', QdiffF);
fprintf('================================\n');

% Time-domain plots.
figure('Color','w','Name','Two GFL theta stability probe');
tiledlayout(3,2,'TileSpacing','compact','Padding','compact');
sgtitle(sprintf('Two GFL no-fault probe, initial \\Delta\\theta = %+g deg', thetaOffsetDeg));

nexttile([1 2]);
plot(t, thetaDiffDeg, 'LineWidth', 1.2);
grid on;
xlabel('Time (s)');
ylabel('\theta_1 - \theta_2 (deg)');
title('Angle difference');

nexttile;
plot(t, store.Q1, 'LineWidth', 1.1); hold on;
plot(t, store.Q2, 'LineWidth', 1.1);
grid on;
xlabel('Time (s)');
ylabel('Q (pu)');
title('Reactive powers');
legend('Q_1','Q_2','Location','best');

nexttile;
plot(t, Qdiff, 'LineWidth', 1.2);
grid on;
xlabel('Time (s)');
ylabel('Q_1 - Q_2 (pu)');
title('Reactive power difference');

nexttile;
plot(t, store.w1, 'LineWidth', 1.1); hold on;
plot(t, store.w2, 'LineWidth', 1.1);
grid on;
xlabel('Time (s)');
ylabel('\omega (pu)');
title('PLL frequency');
legend('\omega_1','\omega_2','Location','best');

nexttile;
plot(t, store.vd1, 'LineWidth', 1.1); hold on;
plot(t, store.vq1, 'LineWidth', 1.1);
grid on;
xlabel('Time (s)');
ylabel('v_{dq}');
title('GFL1 terminal voltage');
legend('v_{d1}','v_{q1}','Location','best');

% Phase portrait: Delta(theta) versus Delta(Q).
figure('Color','w','Name','DeltaTheta_vs_DeltaQ');
scatter(thetaDiffDeg, Qdiff, 25, t, 'filled');
hold on;
plot(thetaDiffDeg, Qdiff, 'k-', 'LineWidth', 1.0);
plot(thetaDiffDeg(1), Qdiff(1), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 9);
plot(thetaDiffDeg(end), Qdiff(end), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 9);
grid on;
xlabel('\Delta\theta = \theta_1 - \theta_2 (deg)');
ylabel('\DeltaQ = Q_1 - Q_2 (pu)');
title('\Delta\theta-\DeltaQ trajectory');
colorbar;
colormap(parula);
legend('Trajectory samples','Trajectory','Initial','Final','Location','best');

% Phase portrait: Delta(theta) versus Delta(omega).
figure('Color','w','Name','DeltaTheta_vs_DeltaOmega');
scatter(thetaDiffDeg, deltaOmega, 25, t, 'filled');
hold on;
plot(thetaDiffDeg, deltaOmega, 'k-', 'LineWidth', 1.0);
plot(thetaDiffDeg(1), deltaOmega(1), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 9);
plot(thetaDiffDeg(end), deltaOmega(end), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 9);
grid on;
xlabel('\Delta\theta = \theta_1 - \theta_2 (deg)');
ylabel('\Delta\omega = \omega_1 - \omega_2 (pu)');
title('\Delta\theta-\Delta\omega trajectory');
colorbar;
colormap(parula);
legend('Trajectory samples','Trajectory','Initial','Final','Location','best');

% I-V trajectory (terminal voltage magnitude vs injected current magnitude).
figure('Color','w','Name','IV_Trajectories');
tiledlayout(1,2,'TileSpacing','compact','Padding','compact');
sgtitle(sprintf('I-V trajectories, initial \Delta\theta = %+g deg', thetaOffsetDeg));

nexttile;
scatter(store.V1, store.I1, 25, t, 'filled');
hold on;
plot(store.V1, store.I1, 'k-', 'LineWidth', 1.0);
plot(store.V1(1), store.I1(1), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 9);
plot(store.V1(end), store.I1(end), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 9);
grid on;
xlabel('|V_1| (pu)');
ylabel('|I_1| (pu)');
title('GFL1 |I| vs |V|');
colorbar;
colormap(parula);

nexttile;
scatter(store.V2, store.I2, 25, t, 'filled');
hold on;
plot(store.V2, store.I2, 'k-', 'LineWidth', 1.0);
plot(store.V2(1), store.I2(1), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 9);
plot(store.V2(end), store.I2(end), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 9);
grid on;
xlabel('|V_2| (pu)');
ylabel('|I_2| (pu)');
title('GFL2 |I| vs |V|');
colorbar;
colormap(parula);

end

function dx = ode_case(t, x, p, caseCfg)
    s = state_to_struct(x); %#ok<NASGU>
    out = solve_two_bus_network_gfl_avg_pi(s, p, caseCfg, false);

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
    % Fixed-point alignment to the pre-perturbation equilibrium.
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
    store.Q1  = zeros(N,1); store.Q2  = zeros(N,1);
    store.I1  = zeros(N,1); store.I2  = zeros(N,1);
    store.V1  = zeros(N,1); store.V2  = zeros(N,1);
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
    store.Q1(k)  = out.Q1;  store.Q2(k)  = out.Q2;
    store.I1(k)   = abs(out.I1); store.I2(k) = abs(out.I2);
    store.V1(k)   = abs(out.V1); store.V2(k) = abs(out.V2);
end

function a = wrap_pi(x)
    a = atan2(sin(x), cos(x));
end

function out = solve_two_bus_network_gfl_avg_pi(x, p, caseCfg, faultActive)
%SOLVE_TWO_BUS_NETWORK_GFL_AVG_PI
% Two-bus resistance-dominated network for the averaged two-GFL case.
%
% The inverters are represented by current states (dq in each inverter's
% own rotating frame). Their injected currents are solved against the
% network admittance.

th1 = x.th1;
th2 = x.th2;

switch caseCfg.type1
    case 'GFL'
        I1dq = x.Id1 + 1j*x.Iq1;
        Iinj1 = I1dq * exp(1j*th1);
    otherwise
        error('solve_two_bus_network_gfl_avg_pi only supports GFL-GFL.');
end

switch caseCfg.type2
    case 'GFL'
        I2dq = x.Id2 + 1j*x.Iq2;
        Iinj2 = I2dq * exp(1j*th2);
    otherwise
        error('solve_two_bus_network_gfl_avg_pi only supports GFL-GFL.');
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

Yfault = 0;
if faultActive
    Yfault = p.Yfault;
end

Ybus = [Y1 + Y12,       -Y12;
             -Y12, Y2 + Y12 + Yfault];

V = Ybus \ [Iinj1; Iinj2];
V1 = V(1);
V2 = V(2);

% dq variables in each inverter's own rotating frame.
v1dq = V1 * exp(-1j*th1);
v2dq = V2 * exp(-1j*th2);
i1dq = Iinj1 * exp(-1j*th1);
i2dq = Iinj2 * exp(-1j*th2);

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

function p = paper_two_gfl_common_avg_pi()
%PAPER_TWO_GFL_COMMON_AVG_PI
% Parameters for a reduced-order two-GFL fault case with averaged current
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
% 0.0 gives the strongest fault-current transient; 1.0 gives perfect
% cancellation and very flat current. A small value keeps the response
% physical but not too aggressive.
p.kff = 0.0000001;

% Simple back-calculation anti-windup gain.
p.kaw = 6.0;

% Optional converter-voltage saturation.
p.vConvMax = 3;

% Grid-following current references.
p.Idref1 = 0.5;
p.Iqref1 = 0.5; % 0.0;
p.Idref2 = 0.5;
p.Iqref2 = -0.5; %0

% Network model.
p.Vref = 1.0;
p.Gload1 = 0.4;
p.Gload2 = 0.6;
p.Z12    = 0.5 + 1j*0.3;
p.Yfault = 1e4;

% Optional overall current ceiling used only as a safety check.
p.iMax = 5.0;

% Simulation window.
p.tEnd = 0.60;
p.tFaultOn = 0.20;
p.tFaultOff = 0.26;



end




