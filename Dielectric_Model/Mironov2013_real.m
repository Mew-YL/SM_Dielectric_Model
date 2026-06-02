function eps_real = Mironov2013_real(vwc, T, vclay)
% Mironov2013_real: Calculate only the real part of the dielectric constant (ε')
% f      - frequency [GHz] (保留接口，但未显式使用)
% vwc    - volumetric water content [m3/m3]
% T      - temperature [°C]
% vclay  - clay fraction [0–1]

    if nargin ~=3
        error('Usage: eps_real = Mironov2013_real(f, vwc, T, vclay)');
    end

    % --- Convert clay fraction to percent (模型经验公式基于 % clay)
    clay_pct = vclay; 
    
    % --- 最大结合水含量 (mvt)
    mvt = 0.02863 + 0.0030673 * clay_pct;

    % --- 折射率和衰减参数 ---
    nd = 1.634 - 0.00539 * clay_pct + 2.75e-5 * clay_pct.^2;
    kd = 0.0395 - 4.083e-4 * clay_pct;

    nb = 8.86 + 0.00321 * T ...
       + (-0.0644 + 7.96e-4 * T) * clay_pct ...
       + (2.97e-4 - 9.6e-6 * T) * clay_pct.^2;

    kb = 0.738 - 0.00903 * T + 8.57e-5 * T.^2 ...
       + (-0.00215 + 1.47e-4 * T) * clay_pct ...
       + (7.36e-5 - 1.03e-6 * T + 1.05e-8 * T.^2) * clay_pct.^2;

    nu = (10.3 - 0.0173 * T) ...
       + (6.5e-4 + 8.82e-5 * T) * clay_pct ...
       + (-6.34e-6 - 6.32e-7 * T) * clay_pct.^2;

    ku = (0.7 - 0.017 * T + 1.78e-4 * T.^2) ...
       + (0.0161 + 7.25e-4 * T) * clay_pct ...
       + (-1.146e-4 - 6.03e-6 * T - 7.87e-9 * T.^2) * clay_pct.^2;

    % --- 按体积含水量分段计算 ---
    if vwc <= mvt
        n = nd + (nb - 1) * vwc;
        k = kd + kb * vwc;
    else
        n = nd + (nb - 1) * mvt + (nu - 1) * (vwc - mvt);
        k = kd + kb * mvt + ku * (vwc - mvt);
    end

    % --- 实部计算 ---
    eps_real = n.^2 - k.^2;
end
