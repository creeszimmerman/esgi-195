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
