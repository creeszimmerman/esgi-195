function theta_sweep_qdiff_stability()
%THETA_SWEEP_QDIFF_STABILITY
% Sweep initial theta1-theta2 from -180 to +180 degrees.
% For each initial angle offset:
%   1) build the pre-perturbation equilibrium,
%   2) apply the angle split,
%   3) simulate with NO short-circuit/fault,
%   4) save initial/final theta1-theta2 and Q1-Q2 to CSV.
%
% The script also generates a quick summary figure at the end.

clear; close all; clc;

p = make_params();
caseCfg.type1 = 'GFL';
caseCfg.type2 = 'GFL';

% Sweep resolution.
deltaVecDeg = -180:5:180;
nCases = numel(deltaVecDeg);

% Output files.
csvFile = fullfile(pwd, 'theta_sweep_results.csv');
figFile = fullfile(pwd, 'theta_sweep_summary.png');

% Solver options.
opts = odeset('RelTol',1e-6,'AbsTol',1e-8,'MaxStep',1e-3);

% Build the nominal equilibrium once.
xEq = build_equilibrium_state(p, caseCfg);

% Preallocate results.
results = table('Size',[nCases 9], ...
    'VariableTypes',{'double','double','double','double','double','double','logical','string','double'}, ...
    'VariableNames',{'theta0_deg','thetaFinal_deg','Qdiff0','QdiffFinal','Q1_0','Q2_0','converged','status','tEnd'});

% Store some basic summary metrics for plots.
theta0_all = nan(nCases,1);
thetaF_all  = nan(nCases,1);
Qdiff0_all  = nan(nCases,1);
QdiffF_all  = nan(nCases,1);
ok_all      = false(nCases,1);

fprintf('Running %d cases from %d to %d deg...\n', nCases, deltaVecDeg(1), deltaVecDeg(end));

for k = 1:nCases
    deltaDeg = deltaVecDeg(k);
    delta = deg2rad(deltaDeg);

    % Start from the nominal equilibrium, then split the angles.
    x0 = xEq;
    x0.th1 = wrap_pi(x0.th1 + delta/2);
    x0.th2 = wrap_pi(x0.th2 - delta/2);

    % Compute the initial Q values directly from the perturbed initial state.
    try
        out0 = solve_two_bus_network_gfl_avg_pi_no_fault(x0, p, caseCfg);
        theta0 = wrap_pi(x0.th1 - x0.th2);
        Qdiff0 = out0.Q1 - out0.Q2;
        Q1_0 = out0.Q1;
        Q2_0 = out0.Q2;
    catch ME
        fprintf('Case %d/%d, delta=%g deg: initial network solve failed (%s)\n', k, nCases, deltaDeg, ME.message);
        theta0 = NaN;
        Qdiff0 = NaN;
        Q1_0 = NaN;
        Q2_0 = NaN;
    end

    % Simulate.
    converged = false;
    status = "ok";
    thetaF = NaN;
    QdiffF = NaN;
    tEndReached = NaN;

    try
        xInit = struct_to_state(x0);
        [tSim, xSim] = ode15s(@(t,x) ode_case_no_fault(t, x, p, caseCfg), [0 p.tEnd], xInit, opts);

        xEnd = state_to_struct(xSim(end,:).');
        outEnd = solve_two_bus_network_gfl_avg_pi_no_fault(xEnd, p, caseCfg);

        thetaF = wrap_pi(xEnd.th1 - xEnd.th2);
        QdiffF = outEnd.Q1 - outEnd.Q2;
        tEndReached = tSim(end);

        % Rough convergence screen.
        converged = isfinite(thetaF) && isfinite(QdiffF) && ...
                    abs(rad2deg(thetaF)) < 5 && ...
                    abs(QdiffF) < 1e-2;

    catch ME
        status = "fail";
        fprintf('Case %d/%d, delta=%g deg: simulation failed (%s)\n', k, nCases, deltaDeg, ME.message);
    end

    theta0_all(k) = rad2deg(theta0);
    thetaF_all(k)  = rad2deg(thetaF);
    Qdiff0_all(k)  = Qdiff0;
    QdiffF_all(k)  = QdiffF;
    ok_all(k)      = converged;

    results.theta0_deg(k)    = theta0_all(k);
    results.thetaFinal_deg(k)= thetaF_all(k);
    results.Qdiff0(k)        = Qdiff0_all(k);
    results.QdiffFinal(k)    = QdiffF_all(k);
    results.Q1_0(k)          = Q1_0;
    results.Q2_0(k)          = Q2_0;
    results.converged(k)     = converged;
    results.status(k)        = status;
    results.tEnd(k)          = tEndReached;

    fprintf('%4d/%4d  theta0=%7.2f deg  thetaF=%8.2f deg  Qdiff0=% .5f  QdiffF=% .5f  %s\n', ...
        k, nCases, theta0_all(k), thetaF_all(k), Qdiff0_all(k), QdiffF_all(k), status);
end

writetable(results, csvFile);
fprintf('\nWrote CSV: %s\n', csvFile);

% Summary figure.
figure('Color','w','Name','Theta sweep summary');
tiledlayout(2,2,'TileSpacing','compact','Padding','compact');

nexttile;
plot(theta0_all, thetaF_all, '.', 'MarkerSize', 10);
grid on;
xlabel('Initial \theta_1 - \theta_2 (deg)');
ylabel('Final \theta_1 - \theta_2 (deg)');
title('Angle map');

nexttile;
plot(theta0_all, Qdiff0_all, '.', 'MarkerSize', 10);
grid on;
xlabel('Initial \theta_1 - \theta_2 (deg)');
ylabel('Initial Q_1 - Q_2');
title('Initial reactive-power difference');

nexttile;
plot(theta0_all, QdiffF_all, '.', 'MarkerSize', 10);
grid on;
xlabel('Initial \theta_1 - \theta_2 (deg)');
ylabel('Final Q_1 - Q_2');
title('Final reactive-power difference');

nexttile;
stairs(theta0_all, double(ok_all), 'LineWidth', 1.2);
grid on;
xlabel('Initial \theta_1 - \theta_2 (deg)');
ylabel('Converged (1=yes)');
title('Rough convergence flag');

saveas(gcf, figFile);
fprintf('Saved figure: %s\n', figFile);

end

function p = make_params()
% Same reduced-order two-GFL setup, with the finite loads restored.
p.f0 = 50;
p.w0 = 2*pi*p.f0;

% PLL parameters.
p.wpll  = 2*pi*15;
p.kppll = p.wpll;
p.kipll = p.wpll^2/4;

% Current-control / output-filter surrogate.
p.Tlf = 0.018;
p.Rf  = 0.03;
p.Xf  = 0.12;

% PI current controller gains.
p.kpc = 5;
p.kic = 1000;

% Partial voltage feedforward coefficient.
p.kff = 0.0000001;

% Anti-windup gain.
p.kaw = 6.0;

% Converter-voltage saturation.
p.vConvMax = 3;

% Grid-following current references.
p.Idref1 = 0.5;
p.Iqref1 = 0.5; %0
p.Idref2 = 0.5;
p.Iqref2 = -0.5; %0

% Network model.
p.Vref = 1.0;
p.Gload1 = 0.4;
p.Gload2 = 0.6;
p.Z12    = 0.5 + 1j*0.3;

% No fault in this experiment.
p.Yfault = 0;

% Optional current ceiling.
p.iMax = 5.0;

% Simulation window.
p.tEnd = 2.0;
end

function xEq = build_equilibrium_state(p, caseCfg)
% Fixed-point alignment to the no-fault equilibrium.
xEq = struct();
xEq.th1  = 0;
xEq.xi1  = 0;
xEq.Id1  = p.Idref1;
xEq.Iq1  = p.Iqref1;
xEq.zId1 = 0;
xEq.zIq1 = 0;
xEq.th2  = 0;
xEq.xi2  = 0;
xEq.Id2  = p.Idref2;
xEq.Iq2  = p.Iqref2;
xEq.zId2 = 0;
xEq.zIq2 = 0;

for iter = 1:40 %#ok<NASGU>
    out = solve_two_bus_network_gfl_avg_pi_no_fault(xEq, p, caseCfg);

    % Align PLL frames to the local bus voltages.
    xEq.th1 = angle(out.V1);
    xEq.th2 = angle(out.V2);

    % Keep the current references at their commanded values.
    xEq.Id1 = p.Idref1;
    xEq.Iq1 = p.Iqref1;
    xEq.Id2 = p.Idref2;
    xEq.Iq2 = p.Iqref2;

    % Steady-state controller bias.
    v1dq = out.V1 * exp(-1j * xEq.th1);
    v2dq = out.V2 * exp(-1j * xEq.th2);
    v1d = real(v1dq);
    v1q = imag(v1dq);
    v2d = real(v2dq);
    v2q = imag(v2dq);

    xEq.zId1 = (1 - p.kff) * v1d;
    xEq.zIq1 = (1 - p.kff) * v1q;
    xEq.zId2 = (1 - p.kff) * v2d;
    xEq.zIq2 = (1 - p.kff) * v2q;

    % PLL integrators so omega starts nominal.
    xEq.xi1 = -p.kppll * v1q;
    xEq.xi2 = -p.kppll * v2q;
end
end

function dx = ode_case_no_fault(~, x, p, caseCfg)
s = state_to_struct(x);
out = solve_two_bus_network_gfl_avg_pi_no_fault(s, p, caseCfg);

% Inverter 1
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

function out = solve_two_bus_network_gfl_avg_pi_no_fault(x, p, caseCfg)
% Same topology as the solver in your folder, but with Yfault = 0.
th1 = x.th1;
th2 = x.th2;

switch caseCfg.type1
    case 'GFL'
        I1dq = x.Id1 + 1j*x.Iq1;
        Iinj1 = I1dq * exp(1j*th1);
    otherwise
        error('Only GFL-GFL is supported.');
end

switch caseCfg.type2
    case 'GFL'
        I2dq = x.Id2 + 1j*x.Iq2;
        Iinj2 = I2dq * exp(1j*th2);
    otherwise
        error('Only GFL-GFL is supported.');
end

Y12 = 1 / p.Z12;
Y1  = p.Gload1;
Y2  = p.Gload2;
Yfault = 0;

Ybus = [Y1 + Y12,       -Y12;
             -Y12, Y2 + Y12 + Yfault];

V = Ybus \ [Iinj1; Iinj2];
V1 = V(1);
V2 = V(2);

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
end

function x = struct_to_state(s)
x = [s.th1; s.xi1; s.Id1; s.Iq1; s.zId1; s.zIq1; ...
     s.th2; s.xi2; s.Id2; s.Iq2; s.zId2; s.zIq2];
end

function a = wrap_pi(x)
a = atan2(sin(x), cos(x));
end
