function dcsr = Mironov2009real(f, vwc, vclay)
% Mironov2019: Calculate real part of soil dielectric constant (ε') 
% f [GHz], vwc [m3/m3], vclay [fraction 0–1]

    if nargin ~= 3
        error('Usage: dcsr = Mironov2099(f, vwc, vclay)');
    end

    f = f * 1e9;           % 转为 Hz
   

    mvt = 0.02863 + 0.0030673 * vclay;
    nd  = 1.634 - 0.00539 * vclay + 2.748e-5 * vclay.^2;
    kd  = 0.03952 - 4.083e-4 * vclay;

    e0b = 79.8 - 0.854 * vclay + 0.00327 * vclay.^2;
    t1b = 1.062e-11 + 3.45e-14 * vclay;
    ob  = 0.3112 + 0.00467 * vclay;

    e0u = 100;
    tu  = 8.5e-12;
    ou  = 0.3631 + 0.01217 * vclay;

    e0 = 8.854187817e-12;

    eb1 = 4.9 + (e0b - 4.9) ./ (1 + (2*pi*f.*t1b).^2);
    eu1 = 4.9 + (e0u - 4.9) ./ (1 + (2*pi*f.*tu).^2);
    eb2 = (2*pi*f.*t1b.*(e0b - 4.9)) ./ (1 + (2*pi*f.*t1b).^2) + ob./(2*pi*f*e0);
    eu2 = (2*pi*f.*tu.*(e0u - 4.9)) ./ (1 + (2*pi*f.*tu).^2) + ou./(2*pi*f*e0);

    nb = sqrt(0.5) * sqrt(sqrt(eb1.^2 + eb2.^2) + eb1);
    kb = sqrt(0.5) * sqrt(sqrt(eb1.^2 + eb2.^2) - eb1);
    nu = sqrt(0.5) * sqrt(sqrt(eu1.^2 + eu2.^2) + eu1);
    ku = sqrt(0.5) * sqrt(sqrt(eu1.^2 + eu2.^2) - eu1);

    if vwc <= mvt
        n = nd + (nb - 1) * vwc;
        k = kd + kb * vwc;
    else
        n = nd + (nb - 1) * mvt + (nu - 1) * (vwc - mvt);
        k = kd + kb * mvt + ku * (vwc - mvt);
    end

    dcsr = n.^2 - k.^2;  % 实部
end
