
function p = paper_gfl_gfm_exact_common()
%PAPER_GFL_GFM_EXACT_COMMON
% Shared parameters for a mixed GFL-GFM reduced-order fault case.
%
% GFL side: PLL + averaged PI current control + explicit current plant.
% GFM side: P-w droop + filtered active power + virtual inductance.
% The network is a two-bus system with an inductive tie line and a
% temporary short-circuit at bus 2.

p.f0 = 50;
p.w0 = 2*pi*p.f0;

% --- GFL (bus 1) ---
p.wpll  = 2*pi*15;
p.kppll = p.wpll;
p.kipll = p.wpll^2/4;

% Explicit averaged current plant.
p.Tlf = 0.02;
p.Rf  = 0.02;
p.Xf  = 0.12;

% Averaged PI current controller.
p.kpc = 4.0;
p.kic = 120.0;
p.kaw = 2.0;

% Partial terminal-voltage feedforward.
p.kff = 0.20;

% GFL current references.
p.Idref1 = 0.50;
p.Iqref1 = 0.00;

% Current-command magnitude limit.
p.vConvMax = 1.8;

% --- GFM (bus 2) ---
p.Vref    = 1.0;
p.mDroop2 = 0.06;
p.wf      = 2*pi*15;
p.Pref2   = 0.20;
p.Zvirt2  = 1j*0.25;

% --- Network ---
p.Gload1 = 0.50;
p.Gload2 = 0.50;
p.Z12    = 1j*0.35;

% Temporary short-circuit at bus 2.
p.Yfault = 1e5;

% Timing for the fault case.
p.tEndFault  = 0.60;
p.tFaultOn   = 0.20;
p.tFaultOff  = 0.26;

% Timing for the Fig. 18-style PLL-freeze case.
p.tEndFig18    = 5.0;
p.tFreezePLL   = 1.0;

% Solver tolerances.
p.relTol = 1e-6;
p.absTol = 1e-8;
p.maxStep = 1e-3;
end
