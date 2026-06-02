function eps_real = mironov2019SOC_real(theta_v, rho_d, T_C, SOM)
% MIRONOV2019_THAWED_REAL
% 计算有机土壤在1.4 GHz的实部介电常数 (thawed状态)
% 基于Mironov et al. 2019 TGRS模型
%
% INPUTS
%   theta_v - 体积含水量 (cm^3/cm^3), 向量或标量
%   rho_d   - 干容重 (g/cm^3), 标量
%   T_C     - 土壤温度 (°C), 标量 (范围: 0...25)
%   SOM     - 土壤有机质含量 (% by weight), 标量 (35..80)
%
% OUTPUT
%   eps_real - 实部介电常数

    % ---- 1) 将体积含水量转换为质量含水量 ----
    % theta_v (cm^3/cm^3) -> mg (g/g)
    rho_w = 1.0;  % 水的密度 g/cm^3
    mg = (rho_w * theta_v) ./ rho_d;
    mg = mg(:);  % 转换为列向量

    % ---- 2) 计算mg1, mg2（仅thawed状态） ----
    [mg1, mg2] = mir19_mg1_mg2_thawed(T_C, SOM);

    % ---- 3) 获取thawed状态参数 ----
    params = mir19_params_thawed(T_C);

    A_no = params.no_1_over_rho_o;       % (no - 1)/rho_o
    A_nb = params.nb_1_over_rho_b;       % (nb - 1)/rho_b
    A_nt = params.nt_1_over_rho_t;       % (nt - 1)/rho_t
    A_nL = params.nL_1_over_rho_L;       % (nl - 1)/rho_l

    K_o  = params.kappa_o_over_rho_o;    % kappa_o/rho_o
    K_b  = params.kappa_b_over_rho_b;    % kappa_b/rho_b
    K_t  = params.kappa_t_over_rho_t;    % kappa_t/rho_t
    K_L  = params.kappa_L_over_rho_L;    % kappa_l/rho_l

    % ---- 4) Reduced CRI: (ns - 1)/rho_d, kappa_s/rho_d，按 mg 分段 ----
    R_n = zeros(size(mg));   % (ns - 1)/rho_d
    R_k = zeros(size(mg));   %  kappa_s/rho_d

    idx1 = (mg <= mg1);
    idx2 = (mg > mg1) & (mg <= mg2);
    idx3 = (mg > mg2);

    % 区间 1: 0 <= mg <= mg1
    R_n(idx1) = A_no + A_nb .* mg(idx1);
    R_k(idx1) = K_o  + K_b  .* mg(idx1);

    % 区间 2: mg1 < mg <= mg2
    R_n(idx2) = A_no + A_nb * mg1 + A_nt .* (mg(idx2) - mg1);
    R_k(idx2) = K_o  + K_b  * mg1 + K_t  .* (mg(idx2) - mg1);

    % 区间 3: mg > mg2
    R_n(idx3) = A_no + A_nb * mg1 + A_nt * (mg2 - mg1) ...
                      + A_nL .* (mg(idx3) - mg2);
    R_k(idx3) = K_o  + K_b  * mg1 + K_t  * (mg2 - mg1) ...
                      + K_L  .* (mg(idx3) - mg2);

    % ---- 5) 从 reduced CRI 还原 ns, kappa ----
    ns    = 1 + rho_d .* R_n;
    kappa =      rho_d .* R_k;

    % ---- 6) 只计算并返回实部介电常数 ----
    eps_real = ns.^2 - kappa.^2;
end

%% ===== 子函数: mg1, mg2 (仅thawed状态) =====
function [mg1, mg2] = mir19_mg1_mg2_thawed(T, SOM)
% 计算thawed状态的mg1和mg2
    % mg1 = 0.118 + 8.695×10^-4 * SOM - 9.6×10^-4 * T
    mg1 = 0.118 + 8.695e-4 * SOM - 9.6e-4 * T;

    % mg2 = 0.382 + 9.208×10^-4 * SOM - 1.91×10^-3 * T
    mg2 = 0.382 + 9.208e-4 * SOM - 1.91e-3 * T;
end

%% ===== 子函数: 参数获取 (仅thawed状态) =====
function params = mir19_params_thawed(T)
% 获取thawed状态的所有参数
    % (no - 1)/rho_o = 0.504 + 8.75×10^-7 * T
    params.no_1_over_rho_o = 0.504 + 8.75e-7 * T;

    % (nb - 1)/rho_b = 3.010 + 0.0328 * T
    params.nb_1_over_rho_b = 3.010 + 0.0328 * T;

    % (nt - 1)/rho_t = 7.572 - 8.33×10^-4 * T
    params.nt_1_over_rho_t = 7.572 - 8.33e-4 * T;

    % (nl - 1)/rho_l = 8.906 - 0.0207 * T
    params.nL_1_over_rho_L = 8.906 - 0.0207 * T;

    % kappa_o / rho_o = 0
    params.kappa_o_over_rho_o = 0.0;

    % kappa_b / rho_b = 1.057 + 2.39×10^-3 * T
    params.kappa_b_over_rho_b = 1.057 + 2.39e-3 * T;

    % kappa_t / rho_t = 1.831 - 0.0252*T
    params.kappa_t_over_rho_t = 1.831 - 0.0252 * T;

    % kappa_l / rho_l = 0.832 - 2.21×10^-2*T + 4.37×10^-4*T^2
    params.kappa_L_over_rho_L = 0.832 - 2.21e-2*T + 4.37e-4*T.^2;
end
