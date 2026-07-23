function p = paper_two_gfm_common()
%PAPER_TWO_GFM_COMMON
% Parameters for a reduced-order two-GFM fault case with P-w droop control.
%
% The model uses two voltage-forming sources behind virtual inductances,
% connected by an inductive tie line, with a temporary short at bus 2.

p.f0 = 50;
p.w0 = 2*pi*p.f0;

% Nominal voltage magnitude of the internal GFM sources.
p.Vref = 1.0;

% Active-power droop (per-unit frequency slope).
% omega_pu = 1 + mDroop*(Pref - Pfiltered)
p.mDroop1 = 0.06;
p.mDroop2 = 0.06;

% Low-pass filter on measured active power.
p.wf = 2*pi*15;

% Symmetric power transfer so the pre-fault operating point is well behaved.
p.Pref1 = 0.20;
p.Pref2 = -0.20;

% Virtual impedance of each GFM output stage (purely inductive for the
% reduced-order study; this avoids resistive loss imbalance).
p.Zvirt1 = 1j*0.25;
p.Zvirt2 = 1j*0.25;

% Inductive connection between the two buses.
p.Z12 = 1j*0.35;

% Temporary short-circuit at bus 2.
p.Yfault = 1e5;

% Simulation window.
p.tEnd = 0.60;
p.tFaultOn = 0.20;
p.tFaultOff = 0.26;

% ODE solver tolerances.
p.relTol = 1e-6;
p.absTol = 1e-8;
p.maxStep = 1e-3;
end
