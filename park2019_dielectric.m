function [eps_real, eps_imag, eps_complex] = park2019_dielectric( ...
    theta, v_clay, v_silt, v_sand, OM, ...
    freq_Hz, T_C, ECe, BD, porosity, w_wp)
%PARK2019_DIELECTRIC  Complex effective dielectric constant (Park 2019 model)
%
% 输入参数
%   theta    : 体积含水量 w (m^3/m^3)
%   v_clay   : 黏粒体积分数 (cm^3/cm^3) 0-1
%   v_silt   : 粉粒体积分数 (cm^3/cm^3) 0-1
%   v_sand   : 砂粒体积分数 (cm^3/cm^3) 0-1
%   OM       : 有机质含量 OM (%)，质量百分数
%   freq_Hz  : 频率 Hz，缺省 1.4e9
%   T_C      : 土壤温度 ℃，缺省 20
%   ECe      : 盐度 ECe (dS/m)，缺省 0
%   BD       : 土壤容重 BD (g/cm^3)，若为空则用 Park 2019 式(3) 计算
%   porosity : 孔隙度 p (m^3/m^3)，若为空则用 Eq.(2) 计算
%   w_wp     : 萎蔫点 w_wp (m^3/m^3)，若为空则用 Eq.(1) 计算
%
% 输出参数
%   eps_real    : 有效介电常数实部 ε'
%   eps_imag    : 有效介电常数虚部 ε'' (>0)
%   eps_complex : 复介电常数 ε = ε' - j ε''
%
% 参考：Park et al., 2017, 2019; Toth et al., 2015

    if nargin < 6 || isempty(freq_Hz), freq_Hz = 1.4e9; end
    if nargin < 7 || isempty(T_C),     T_C     = 20.0;  end
    if nargin < 8 || isempty(ECe),     ECe     = 0.0;   end

    % -------- 1) 归一化纹理 --------
    v_sum  = v_clay + v_silt + v_sand;
    if v_sum <= 0
        error('v_clay + v_silt + v_sand 必须 > 0');
    end
    v_clay = v_clay ./ v_sum;
    v_silt = v_silt ./ v_sum;
    v_sand = v_sand ./ v_sum;

    % -------- 2) ECe 转换为盐度 S (‰) --------
    S_ppt = ECe_to_salinity(ECe);

    % -------- 3) 体积质量 BD (Park 2019 Eq.(3)) --------
    if nargin < 9 || isempty(BD)
        BD = bulk_density_park2019(OM);   % g/cm3
    end

    % -------- 4) 孔隙度 p (Park 2019 Eq.(2), 引用 Tóth 2015) --------
    if nargin < 10 || isempty(porosity)
        porosity = porosity_park2019(v_clay, v_silt, BD, OM);
        % 物理上限制在 [0.05, 0.9] 之间
        porosity = max(0.05, min(0.9, porosity));
    end

    % -------- 5) 萎蔫点 w_wp (Park 2019 Eq.(1)) --------
    if nargin < 11 || isempty(w_wp)
        % 注意：这里 v_clay 为体积分数 (cm3/cm3)，OM 为 %
        w_wp = 0.02982 * v_clay + 0.089 + 0.00786 * OM;
    end

    % -------- 6) 各相复介电常数 --------
    eps_air   = complex(1.0, 0.0);                      % 空气
    eps_solid = eps_soil_mineral(v_sand, v_silt, v_clay);  % 矿物
    eps_free  = debye_water_free(freq_Hz, T_C, S_ppt);     % 自由水
    eps_bound = debye_water_bound(freq_Hz, v_clay);        % 束缚水

    w   = theta;      % 体积含水量
    p   = porosity;
    wwp = w_wp;

    % -------- 7) 三个含水区间 (Park 2017 式 (10),(16),(20)) --------
    if w <= wwp
        % 干燥区：只有束缚水
        % ε_eff = (1-p) ε_solid + w ε_bound + (p - w) ε_air
        eps_complex = (1 - p) * eps_solid + w * eps_bound + (p - w) * eps_air;

    elseif w <= p && (p - wwp) > 1e-12
        % 过渡区：束缚水 + 自由水
        % v_bound = (p - w)/(p - w_wp)
        % v_free  = (w - w_wp)/(p - w_wp)
        v_bound   = (p - w)   / (p - wwp);
        v_free    = (w - wwp) / (p - wwp);
        eps_water = v_bound * eps_bound + v_free * eps_free;

        % ε_eff = (1-p) ε_solid + w ε_water + (p - w) ε_air
        eps_complex = (1 - p) * eps_solid + w * eps_water + (p - w) * eps_air;

    else
        % 饱和以上区：ε_eff = (1 - w) ε_solid + w ε_free
        eps_complex = (1 - w) * eps_solid + w * eps_free;
    end

    eps_real = real(eps_complex);
    eps_imag = -imag(eps_complex);  % 约定 ε = ε' - j ε''，输出正的 ε''

end


%% --- 子函数：孔隙度 (Park 2019 Eq.(2)) -----------------------------------
function p = porosity_park2019(v_clay, v_silt, BD, OM)
% v_clay, v_silt : 体积分数 (cm3/cm3)
% BD             : g/cm3
% OM             : 有机质含量 (%)
% OC             : 有机碳含量 (%), OC = OM / 1.724

    OC = OM / 1.724;

    p = 0.6819 ...
        - 0.06480 ./ (OC + 1.0) ...
        - 0.11900 .* BD.^2 ...
        - 0.02668 ...
        + 0.1489  .* v_clay ...
        + 0.08031 .* v_silt ...
        + 0.02321 ./ (OC + 1.0) .* BD.^2 ...
        + 0.01908 .* BD.^2 ...
        - 0.11090 .* v_clay ...
        - 0.2315  .* v_silt .* v_clay ...
        - 0.01197 .* v_silt .* BD.^2 ...
        - 0.01068 .* v_clay .* BD.^2;
end


%% --- 子函数：体积质量 BD (Park 2019 Eq.(3)) ------------------------------
function BD = bulk_density_park2019(OM)
% OM : 有机质含量 (%)
    BD = -0.039 * OM + 1.2301;  % g/cm3
end


%% --- 子函数：自由水 Debye 模型 (Eq.(25)-(31)) ----------------------------
function eps = debye_water_free(freq_Hz, T_C, S)
% freq_Hz : 频率 (Hz)
% T_C     : 温度 (°C)
% S       : 盐度 (‰)

    T = T_C;
    eps_min = 4.9;  % Eq.(27)

    % Eq.(28) ε_max0(T)
    eps_max_0 = 88.045 - 0.4147*T + 6.295e-4*T.^2 + 1.075e-5*T.^3;

    % Eq.(29) a_ST(S,T)
    a_ST = 1.0 + 1.613e-3*S.*T - 3.656e-3*S + 3.21e-5*S.^2 - 4.232e-7*S.^3;
    eps_max = eps_max_0 .* a_ST;

    % Eq.(31) b_ST(S,T)
    b_ST = 1.0 + 2.282e-5*S.*T - 7.638e-4*S - 7.760e-6*S.^2 + 1.105e-8*S.^3;

    % Eq.(30) τ_free(T,S)，注意最后除以 2π
    tau0 = 1.1109e-10 ...
         - 3.824e-12*T ...
         + 6.938e-14*T.^2 ...
         - 5.096e-16*T.^3;
    tau = (tau0 .* b_ST) ./ (2*pi);

    omega = 2*pi*freq_Hz;
    x = omega .* tau;

    real_part = eps_min + (eps_max - eps_min) ./ (1 + x.^2);
    imag_part = x .* (eps_max - eps_min) ./ (1 + x.^2);

    eps = complex(real_part, -imag_part);  % ε = ε' - j ε''
end


%% --- 子函数：束缚水 Debye 模型 (Eq.(25),(26),(33)) -----------------------
function eps = debye_water_bound(freq_Hz, v_clay)
% v_clay : 黏粒体积分数 (cm3/cm3)

    eps_min = 4.9;
    eps_max = 36.0 .* v_clay + 44.0;  % Eq.(32)

    tau = 1e-11;  % s, Eq.(33)
    omega = 2*pi*freq_Hz;
    x = omega .* tau;

    real_part = eps_min + (eps_max - eps_min) ./ (1 + x.^2);
    imag_part = x .* (eps_max - eps_min) ./ (1 + x.^2);

    eps = complex(real_part, -imag_part);
end


%% --- 子函数：矿物相介电常数 (Table 6) ------------------------------------
function eps_solid = eps_soil_mineral(v_sand, v_silt, v_clay)
% 使用 Table 6:
%   sand: ε' = 3,  ε'' = 0.078
%   silt: ε' = 5,  ε'' = 0.078
%   clay: ε' = 5,  ε'' = 0.078

    s = v_sand + v_silt + v_clay;
    if s <= 0
        error('纹理体积分数之和必须 > 0');
    end
    v_sand = v_sand ./ s;
    v_silt = v_silt ./ s;
    v_clay = v_clay ./ s;

    eps_sand = complex(3.0, 0.078);
    eps_silt = complex(5.0, 0.078);
    eps_clay = complex(5.0, 0.078);

    eps_solid = v_sand .* eps_sand + v_silt .* eps_silt + v_clay .* eps_clay;
end


%% --- 子函数：ECe 转换为盐度 ----------------------------------------------
function S_ppt = ECe_to_salinity(ECe_dSm)
%ECe_TO_SALINITY  把 ECe (dS/m) 近似转换成盐度 S (‰, g/kg)
%
% 输入
%   ECe_dSm : 饱和土浆浸提液电导率 ECe，单位 dS/m
%             ECe < 5 用 640, ECe >= 5 用 800
%
% 输出
%   S_ppt   : 盐度 S，单位 ‰ (g/kg)，可直接作为 Park 模型的 S

    % 低盐用 640，高盐用 800
    k = 640 * ones(size(ECe_dSm));
    k(ECe_dSm >= 5) = 800;

    % 近似 TDS = k * ECe (mg/L)
    TDS_mgL = k .* ECe_dSm;

    % mg/L → g/L → g/kg ≈ ‰
    S_ppt = TDS_mgL ./ 1000;
end
