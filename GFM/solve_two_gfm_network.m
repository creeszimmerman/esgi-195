function out = solve_two_gfm_network(x, p, caseCfg, faultActive)
%SOLVE_TWO_GFM_NETWORK
% Two-GFM nodal network with virtual impedances and an inductive tie line.
%
% State vector x must provide x.th1 and x.th2 (radians).
% The inverters are modelled as voltage sources Vref*exp(j*theta) behind
% virtual impedances Zvirt1 and Zvirt2.

th1 = x.th1;
th2 = x.th2;

% Internal source voltages.
E1 = p.Vref * exp(1j*th1);
E2 = p.Vref * exp(1j*th2);

Yv1 = 1 / p.Zvirt1;
Yv2 = 1 / p.Zvirt2;
Y12 = 1 / p.Z12;

Yfault = 0;
if faultActive
    Yfault = p.Yfault;
end

% Nodal equations:
%   (Yv1 + Y12) V1 - Y12 V2 = Yv1 * E1
%   -Y12 V1 + (Yv2 + Y12 + Yfault) V2 = Yv2 * E2
Ybus = [Yv1 + Y12,          -Y12;
               -Y12,  Yv2 + Y12 + Yfault];
Ieq = [Yv1 * E1; Yv2 * E2];

V = Ybus \ Ieq;
V1 = V(1);
V2 = V(2);

% Source currents injected into the network through the virtual impedances.
I1 = Yv1 * (E1 - V1);
I2 = Yv2 * (E2 - V2);

% dq-frame quantities in each inverter's own rotating frame.
v1dq = V1 * exp(-1j*th1);
v2dq = V2 * exp(-1j*th2);
i1dq = I1 * exp(-1j*th1);
i2dq = I2 * exp(-1j*th2);

out.V1 = V1;
out.V2 = V2;
out.I1 = I1;
out.I2 = I2;

out.vd1 = real(v1dq);
out.vq1 = imag(v1dq);
out.vd2 = real(v2dq);
out.vq2 = imag(v2dq);
out.id1 = real(i1dq);
out.iq1 = imag(i1dq);
out.id2 = real(i2dq);
out.iq2 = imag(i2dq);

% Active and reactive power at each inverter terminal (PCC side).
out.P1 = real(V1 * conj(I1));
out.Q1 = imag(V1 * conj(I1));
out.P2 = real(V2 * conj(I2));
out.Q2 = imag(V2 * conj(I2));

out.I1mag = abs(I1);
out.I2mag = abs(I2);
out.V1mag = abs(V1);
out.V2mag = abs(V2);
end
