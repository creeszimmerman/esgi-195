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
p.Iqref1 = 0.0;
p.Idref2 = 0.5;
p.Iqref2 = 0.0;

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
