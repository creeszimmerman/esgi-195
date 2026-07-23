
function out = solve_two_bus_network_gfl_gfm_exact(x, p, caseCfg, faultActive)
%SOLVE_TWO_BUS_NETWORK_GFL_GFM_EXACT
% Two-bus mixed system: GFL current source at bus 1 and GFM voltage source
% behind virtual impedance at bus 2.
%
% Bus 1 (GFL) injects the commanded current in its own dq frame.
% Bus 2 (GFM) is represented as an internal voltage source behind Zvirt2.

th1 = x.th1;
th2 = x.th2;

% Grid-following inverter at bus 1.
I1dq = x.Id1 + 1j*x.Iq1;
I1   = I1dq * exp(1j*th1);

% Grid-forming inverter at bus 2.
E2   = p.Vref * exp(1j*th2);
Yv2  = 1 / p.Zvirt2;
Y12  = 1 / p.Z12;

Yfault = 0;
if faultActive
    Yfault = p.Yfault;
end

% Nodal admittance system.
%   (Gload1 + Y12) V1 - Y12 V2 = I1
%   -Y12 V1 + (Gload2 + Y12 + Yv2 + Yfault) V2 = Yv2*E2
Ybus = [p.Gload1 + Y12,                 -Y12;
                   -Y12, p.Gload2 + Y12 + Yv2 + Yfault];
RHS = [I1; Yv2 * E2];

V = Ybus \ RHS;
V1 = V(1);
V2 = V(2);

% Current through the GFM virtual impedance (source current injected).
I2 = Yv2 * (E2 - V2);

% dq variables in each inverter's own frame.
v1dq = V1 * exp(-1j*th1);
v2dq = V2 * exp(-1j*th2);
i1dq = I1 * exp(-1j*th1);
i2dq = I2 * exp(-1j*th2);

out.V1 = V1; out.V2 = V2;
out.I1 = I1; out.I2 = I2;
out.vd1 = real(v1dq); out.vq1 = imag(v1dq);
out.vd2 = real(v2dq); out.vq2 = imag(v2dq);
out.id1 = real(i1dq); out.iq1 = imag(i1dq);
out.id2 = real(i2dq); out.iq2 = imag(i2dq);
out.P1 = real(V1 * conj(I1));
out.Q1 = imag(V1 * conj(I1));
out.P2 = real(V2 * conj(I2));
out.Q2 = imag(V2 * conj(I2));
out.V1mag = abs(V1);
out.V2mag = abs(V2);
out.I1mag = abs(I1);
out.I2mag = abs(I2);
end
