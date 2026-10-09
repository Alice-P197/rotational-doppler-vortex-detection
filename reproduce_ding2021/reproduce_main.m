%% 复现《涡旋光任意入射条件下的旋转物体转速探测》
%  红外与激光工程 2021, 50(9): 20210451  丁友, 丁源圣, 邱松, 刘通, 任元
%
%  复现内容:
%    [1] 公式(10) 任意入射条件下旋转多普勒频移分布 Δf(θ)
%    [2] 文献图3: γ=30°, φ=30°, d=1 mm, r=4 mm, l=±18, f=50 Hz
%    [3] 公式(11)(12) 特征点提取旋转频率方法验证
%    [4] 频谱仿真: 垂直入射单峰 vs 任意入射展宽
%
%  输出图: results/fig3_theory_curve.png
%          results/fig_special_points.png
%          results/fig_spectrum_aligned.png
%          results/fig_spectrum_general.png
%          results/fig_cases_compare.png

clear; clc; close all;

%% ============ 基础参数（文献取值） ============
l     = 18;        % 拓扑荷数（叠加态 ±l）
f_rot = 50;        % 物体旋转频率 (Hz)
Omega = 2*pi*f_rot;% 角速度 (rad/s)
r     = 4e-3;      % 涡旋光半径 (m)
d     = 1e-3;      % 光斑中心横向位移 (m)
gamma = 30*pi/180; % 倾斜角
phi   = 30*pi/180; % 光斑中心-旋转中心连线与长轴夹角

base = l*Omega/pi; % 对准入射频移 lΩ/π = 1800 Hz

fprintf('===== 复现开始 =====\n');
fprintf('基准频移 lΩ/π = %.2f Hz\n', base);

%% ============ [1][2] 复现公式(10) 与文献图3 ============
theta = linspace(0, 2*pi, 2001);

% 公式(10): Δf(θ) = (lΩ/π)[ sin²θ/cosγ + cosγ·cos²θ
%                              + (d/r)(sinθ·cosφ − cosθ·sinφ·cosγ) ]
Df = base .* ( sin(theta).^2 / cos(gamma) ...
             + cos(gamma) .* cos(theta).^2 ...
             + (d/r) .* (sin(theta).*cos(phi) - cos(theta).*sin(phi).*cos(gamma)) );

fprintf('频移范围: [%.2f, %.2f] Hz\n', min(Df), max(Df));

% --- 图3 复现 ---
fig = figure('Position', [100 100 720 460], 'Color', 'w');
plot(theta, Df, 'b-', 'LineWidth', 1.8);
xlim([0 2*pi]); xticks(0:pi/2:2*pi);
xticklabels({'0','π/2','π','3π/2','2π'});
grid on; box on;
xlabel('\theta / rad'); ylabel('Rotational frequency / Hz');
title(sprintf('\\Deltaf(\\theta),  \\gamma=30°, \\phi=30°, d=1mm, r=4mm, l=±18, f=%.0fHz', f_rot));
saveas(fig, fullfile('results', 'fig3_theory_curve.png'));
close(fig);

%% ============ [3] 特征点提取方法验证 (公式11,12) ============
% 特殊角度 θ = 0, π/2, π, 3π/2
Df_0    = base*(cos(gamma) - (d/r)*sin(phi)*cos(gamma));
Df_pi2  = base*(1/cos(gamma) + (d/r)*cos(phi));
Df_pi   = base*(cos(gamma) + (d/r)*sin(phi)*cos(gamma));
Df_3pi2 = base*(1/cos(gamma) - (d/r)*cos(phi));

fprintf('\n特征点频移:\n');
fprintf('  θ=0:     %.2f Hz\n', Df_0);
fprintf('  θ=π/2:   %.2f Hz\n', Df_pi2);
fprintf('  θ=π:     %.2f Hz\n', Df_pi);
fprintf('  θ=3π/2:  %.2f Hz\n', Df_3pi2);

% 公式(12) 提取 f_mod（文献: f_mod = lΩ/π）
% 正确形式(几何平均, 对偶点 d 项精确抵消):
%   f_mod = sqrt( (Δf_0+Δf_π)/2  ×  (Δf_π/2+Δf_3π/2)/2 ) = lΩ/π
f_mod = sqrt( ((Df_0 + Df_pi)/2) * ((Df_pi2 + Df_3pi2)/2) );
f_est = f_mod / (2*l);
fprintf('  f_mod = %.4f Hz, 提取转速 f = %.4f Hz (真值 %.2f Hz, 误差 %.4f%%)\n', ...
    f_mod, f_est, f_rot, abs(f_est-f_rot)/f_rot*100);

% --- 特征点图 ---
fig = figure('Position', [100 100 720 460], 'Color', 'w');
plot(theta, Df, 'b-', 'LineWidth', 1.6); hold on;
sp_theta = [0 pi/2 pi 3*pi/2];
sp_val   = [Df_0 Df_pi2 Df_pi Df_3pi2];
plot(sp_theta, sp_val, 'ro', 'MarkerSize', 8, 'MarkerFaceColor', 'r');
yline(base, 'k--', 'LineWidth', 1.2);
text(0.1, base*1.02, 'l\Omega/\pi', 'FontSize', 11);
xlim([0 2*pi]); xticks(0:pi/2:2*pi);
xticklabels({'0','π/2','π','3π/2','2π'});
grid on; box on;
xlabel('\theta / rad'); ylabel('Rotational frequency / Hz');
title('特征点提取: 红点 = θ=0,π/2,π,3π/2 处频移');
legend('Δf(θ)', '特征点', '基准 lΩ/π', 'Location', 'southeast');
saveas(fig, fullfile('results', 'fig_special_points.png'));
close(fig);

%% ============ [4] 频谱仿真: 对准 vs 任意入射 ============
% 散射点模型: 光斑上 Np 个散射点，θ 均匀分布，
% 每点反射信号频率 = 该点旋转多普勒频移 Δf(θ_k)
% 总信号: s(t) = Σ_k exp(i·2π·Δf(θ_k)·t) → 频谱 = 各频移贡献叠加

Np   = 400;        % 散射点数
T    = 2.0;        % 采样时长 s
Fs   = 50e3;       % 采样率
tvec = (0:round(T*Fs)-1) / Fs;

% --- (a) 对准入射: γ=d=φ=0, Δf = lΩ/π 常数 → 单峰 ---
s_aligned = exp(1i*2*pi*base*tvec);  % 单频

% --- (b) 任意入射: Δf(θ) 分布 → 展宽 ---
th_k = (0:Np-1)/Np * 2*pi;
Df_k = base .* ( sin(th_k).^2 / cos(gamma) ...
               + cos(gamma) .* cos(th_k).^2 ...
               + (d/r) .* (sin(th_k).*cos(phi) - cos(th_k).*sin(phi).*cos(gamma)) );
s_general = sum(exp(1i*2*pi*Df_k(:)*tvec), 1);

% 频谱
f_ax = linspace(0, Fs/2, numel(tvec)/2+1);
fft2amp = @(s) abs(fft(s(:))) / numel(s);

P_aligned = fft2amp(s_aligned);
P_aligned = P_aligned(1:numel(f_ax));
P_general = fft2amp(s_general);
P_general = P_general(1:numel(f_ax));

% 理论频移区间
Df_min = min(Df_k); Df_max = max(Df_k);
fprintf('\n任意入射理论频移区间: [%.2f, %.2f] Hz\n', Df_min, Df_max);

% --- 图: 对准入射频谱 ---
fig = figure('Position', [100 100 720 420], 'Color', 'w');
plot(f_ax, 20*log10(P_aligned/max(P_aligned)+eps), 'b-', 'LineWidth', 1.4);
xlim([base-300 base+300]); ylim([-60 2]);
grid on; box on;
xlabel('Frequency / Hz'); ylabel('Relative amplitude / dB');
title(sprintf('对准入射 (\\gamma=d=\\phi=0): 单峰 f = l\\Omega/\\pi = %.0f Hz', base));
saveas(fig, fullfile('results', 'fig_spectrum_aligned.png'));
close(fig);

% --- 图: 任意入射频谱 ---
fig = figure('Position', [100 100 720 420], 'Color', 'w');
plot(f_ax, 20*log10(P_general/max(P_general)+eps), 'b-', 'LineWidth', 1.4); hold on;
xline(Df_min, 'r--', 'LineWidth', 1.4);
xline(Df_max, 'r--', 'LineWidth', 1.4);
xlim([base-900 base+900]); ylim([-60 2]);
grid on; box on;
xlabel('Frequency / Hz'); ylabel('Relative amplitude / dB');
title(sprintf('任意入射 (\\gamma=30°, \\phi=30°, d=1mm): 频谱展宽 [%.0f, %.0f] Hz', Df_min, Df_max));
legend('频谱', '理论频移边界', 'Location', 'northeast');
saveas(fig, fullfile('results', 'fig_spectrum_general.png'));
close(fig);

%% ============ [5] 文献实验四种工况频移曲线对比 ============
% 原文实验工况（电机设定 53 Hz, l=±18）:
%   工况1(对准): γ≈0,  d=0,   φ=0    实验频移 1906.74 Hz -> 52.96 Hz (<0.1%)
%   工况2:       γ=20°, d=0,   φ=0    实验峰 f_max=2012, f_min=1801
%   工况3:       γ=20°, d=1mm, φ=0    特征: Δf0=Δfπ=1804, f_max=2514, f_min=1538
%   工况4:       γ=20°, d=1mm, φ=30°  特征: 2443/2018/1591/1640
f_rot_exp = 53;                 % 实验电机设定频率 (Hz)
Omega_exp = 2*pi*f_rot_exp;
base_exp  = l*Omega_exp/pi;     % 叠加态对准频移 lΩ/π = 1908 Hz

% 原文实验特征频移（从图5/图6读出）
exp_cases = {
    struct('tag','工况2: γ=20°,d=0,φ=0',    'g',20,'d',0,   'p',0,  'spec',[2012 1801], 'fmod',1903.58, 'fest',52.88, 'err',0.23, 'feat',[1801 2012 1801 2012]);
    struct('tag','工况3: γ=20°,d=1mm,φ=0',  'g',20,'d',1e-3,'p',0,  'spec',[2514 1538], 'fmod',1911.78, 'fest',53.11, 'err',0.21, 'feat',[1804 2514 1804 1538]);
    struct('tag','工况4: γ=20°,d=1mm,φ=30°','g',20,'d',1e-3,'p',30, 'spec',[2443 2018 1591 1640], 'fmod',1919.35, 'fest',53.32, 'err',0.6, 'feat',[2018 2443 1591 1640]);
};

fprintf('\n===== [5] 原文实验四种工况理论验证 =====\n');
fprintf('叠加态对准频移 lΩ/π = %.2f Hz (Ω=2π×53)\n\n', base_exp);

fig = figure('Position', [100 100 860 620], 'Color', 'w');
cols = lines(4);
hold on;
% 工况1 对准
plot(theta, base_exp*ones(size(theta)), '-', 'Color', cols(1,:), 'LineWidth', 1.8, 'DisplayName', '对准 γ≈0');
% 工况2-4 理论曲线 + 实验特征点
for k = 1:numel(exp_cases)
    c = exp_cases{k};
    g_ = c.g*pi/180; d_ = c.d; p_ = c.p*pi/180;
    Df_k2 = base_exp .* ( sin(theta).^2/cos(g_) + cos(g_).*cos(theta).^2 ...
                        + (d_/r).*(sin(theta).*cos(p_) - cos(theta).*sin(p_).*cos(g_)) );
    plot(theta, Df_k2, '-', 'Color', cols(k+1,:), 'LineWidth', 1.6, 'DisplayName', sprintf('\\gamma=%d°, d=%dmm, \\phi=%d°', c.g, c.d*1000, c.p));
    % 理论特征点值
    Df_th = base_exp * [ (cos(g_) - (d_/r)*sin(p_)*cos(g_)), ...
                         (1/cos(g_) + (d_/r)*cos(p_)), ...
                         (cos(g_) + (d_/r)*sin(p_)*cos(g_)), ...
                         (1/cos(g_) - (d_/r)*cos(p_)) ];
    % 实验特征点（红点）
    th_feat = [0 pi/2 pi 3*pi/2];
    for j = 1:numel(c.feat)
        plot(th_feat(j), c.feat(j), 'o', 'Color', cols(k+1,:), 'MarkerSize', 7, 'MarkerFaceColor', cols(k+1,:));
    end
    % 打印理论特征值与实验值对照
    fprintf('%s\n', c.tag);
    fprintf('  理论特征值: θ=0:%.1f  θ=π/2:%.1f  θ=π:%.1f  θ=3π/2:%.1f Hz\n', Df_th(1), Df_th(2), Df_th(3), Df_th(4));
    fprintf('  实验特征值: %s Hz\n', num2str(c.feat, '%.0f '));
    % 复算 f_mod 与转速误差（公式12 几何平均）
    f_mod_c = sqrt( ((c.feat(1)+c.feat(3))/2) * ((c.feat(2)+c.feat(4))/2) );
    f_est_c = f_mod_c/(2*l);
    err_c   = abs(f_est_c - f_rot_exp)/f_rot_exp*100;
    fprintf('  复算 f_mod=%.2f Hz (原文 %.2f), f=%.2f Hz (原文 %.2f), 误差 %.2f%% (原文 %.2f%%)\n\n', ...
        f_mod_c, c.fmod, f_est_c, c.fest, err_c, c.err);
end
xlim([0 2*pi]); xticks(0:pi/2:2*pi);
xticklabels({'0','π/2','π','3π/2','2π'});
grid on; box on;
xlabel('\theta / rad'); ylabel('Rotational frequency / Hz');
title('文献实验四种工况理论频移曲线与实验特征点 (l=±18, f=53Hz)');
legend('Location', 'southeast');
saveas(fig, fullfile('results', 'fig_cases_compare.png'));
close(fig);

%% ============ [6] 对准工况实验验证 ============
% 原文图5: 对准入射实验频移 1906.74 Hz -> f=52.96 Hz (误差<0.1%)
f_detect = 1906.74;
f_est_aligned = f_detect/(2*l);
fprintf('\n===== [6] 对准工况 (原文图5) =====\n');
fprintf('实验频移 %.2f Hz -> f=%.2f Hz (设定 53 Hz, 误差 %.2f%%)\n\n', ...
    f_detect, f_est_aligned, abs(f_est_aligned-53)/53*100);

fprintf('\n===== 复现完成，图片已保存到 results/ =====\n');
