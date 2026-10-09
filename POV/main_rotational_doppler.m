%% 旋转多普勒效应：不同阶数花瓣模光束测量旋转物体速度
%
% 物理原理:
%   花瓣模光束由两个相反 OAM 阶数 (+l, -l) 叠加而成
%   照射旋转物体时，反射光产生频移:
%       Δf_+ = +l·Ω / (2π)
%       Δf_- = -l·Ω / (2π)
%   两分量拍频: f_beat = l·Ω / π
%
%   测量拍频即可反推角速度: Ω_meas = π·f_beat / l
%
% 本脚本仿真不同阶数 l=1~5 的花瓣模光束测量同一旋转目标，
% 对比测量精度与拍频信号质量。
%
% 目录结构:
%   main_rotational_doppler.m      <- 本文件
%   functions/makegrid.m            <- 坐标网格
%   functions/petal_mode_beam.m     <- 花瓣模光束生成
%   functions/rotating_petal_signal.m <- 时域旋转信号仿真
%   functions/extract_beat_freq.m  <- FFT 提取拍频

clear; clc; close all;

% 使用相对于脚本位置的路径，保证从任何目录运行都能找到函数
scriptDir = fileparts(mfilename('fullpath'));
addpath(fullfile(scriptDir, 'functions'));

%% ======================== 参数设置 ========================
% 空间参数
N = 256;               % 采样点数
L = 5e-3;              % 视场大小 (m)
w0 = 0.8e-3;           % 花瓣模束腰半径 (m)

% 旋转目标参数
Omega_true = 2 * pi * 50;   % 真实角速度 (rad/s) → 50 Hz = 3000 rpm
f_true = Omega_true / (2*pi);  % 真实旋转频率 (Hz)

% 仿真时间参数
T_total = 0.2;         % 总仿真时间 (s)
dt = 1e-5;             % 时间步长 (s) → 采样率 100 kHz

% 不同阶数 l
l_list = [1, 2, 3, 4, 5];
num_l = length(l_list);

fprintf('===== 旋转多普勒仿真开始 =====\n');
fprintf('真实旋转频率: %.2f Hz (%.2f rad/s)\n', f_true, Omega_true);
fprintf('仿真时长: %.2f s, 采样率: %.0f Hz\n', T_total, 1/dt);

%% ======================== 生成网格 ========================
[X, Y, R, TH] = makegrid(N, L);

% imagesc 坐标需为向量（取 meshgrid 的第一行/列）
x_ax = X(1,:) * 1e3;   % mm
y_ax = Y(:,1) * 1e3;   % mm

%% ======================== 逐阶数仿真 ========================
results = struct();

for idx = 1:num_l
    l = l_list(idx);
    fprintf('\n--- 阶数 l = %d (花瓣数 = %d) ---\n', l, 2*l);

    % 理论拍频
    f_beat_theory = l * Omega_true / pi;   % Hz
    fprintf('理论拍频: %.2f Hz\n', f_beat_theory);

    % 生成时域信号
    [t, I_sig] = rotating_petal_signal(R, TH, l, w0, Omega_true, T_total, dt, 42);

    % FFT 提取拍频
    f_min = max(1, f_beat_theory * 0.5);
    f_max = f_beat_theory * 1.5;
    [f_axis, P1, f_peak, I_peak] = extract_beat_freq(t, I_sig, f_min, f_max);

    % 反推角速度
    Omega_meas = pi * f_peak / l;   % rad/s
    f_meas = Omega_meas / (2*pi);   % Hz
    error_pct = abs(f_meas - f_true) / f_true * 100;

    fprintf('测量拍频: %.2f Hz\n', f_peak);
    fprintf('反推旋转频率: %.2f Hz (误差 %.2f%%)\n', f_meas, error_pct);

    % 保存结果
    results(idx).l = l;
    results(idx).petals = 2*l;
    results(idx).f_beat_theory = f_beat_theory;
    results(idx).f_peak = f_peak;
    results(idx).f_meas = f_meas;
    results(idx).error_pct = error_pct;
    results(idx).t = t;
    results(idx).I_sig = I_sig;
    results(idx).f_axis = f_axis;
    results(idx).P1 = P1;
    results(idx).f_min = f_min;
    results(idx).f_max = f_max;
end

%% ======================== 画图 ========================

% --- 图1: 不同阶数花瓣模光束强度分布 ---
figure('Name', '花瓣模光束强度分布', 'Position', [100, 100, 1200, 500]);
for idx = 1:num_l
    l = l_list(idx);
    [~, I_petal] = petal_mode_beam(R, TH, l, w0, 0);

    subplot(1, num_l, idx);
    imagesc(x_ax, y_ax, I_petal);
    axis square; axis tight;
    colormap(gca, 'hot');
    xlabel('x (mm)'); ylabel('y (mm)');
    title(sprintf('l = %d, %d petals', l, 2*l));
end
if exist('sgtitle', 'file')
    sgtitle('不同阶数花瓣模光束强度分布');
end

% --- 图2: 时域信号对比 (取前 5 ms) ---
figure('Name', '时域强度信号', 'Position', [100, 100, 1200, 600]);
t_plot_max = 5e-3;  % 只画前 5 ms
for idx = 1:num_l
    subplot(num_l, 1, idx);
    mask = results(idx).t <= t_plot_max;
    plot(results(idx).t(mask)*1e3, results(idx).I_sig(mask), 'LineWidth', 1);
    ylabel(sprintf('l=%d', results(idx).l));
    xlim([0, t_plot_max*1e3]);
    grid on;
    if idx == num_l
        xlabel('时间 (ms)');
    end
end
if exist('sgtitle', 'file')
    sgtitle('不同阶数花瓣模反射光强时域信号（前 5 ms）');
end

% --- 图3: FFT 频谱对比 ---
figure('Name', '拍频频谱', 'Position', [100, 100, 1200, 600]);
for idx = 1:num_l
    subplot(num_l, 1, idx);
    res = results(idx);
    mask = (res.f_axis >= res.f_min*0.8) & (res.f_axis <= res.f_max*1.2);
    plot(res.f_axis(mask), res.P1(mask), 'LineWidth', 1.5);
    hold on;
    if exist('xline', 'file')
        xline(res.f_beat_theory, 'r--', 'LineWidth', 1.5, 'Label', '理论值');
    else
        % 旧版本 MATLAB 用 line 绘制理论值竖线
        yl = ylim;
        hold on;
        line([res.f_beat_theory res.f_beat_theory], yl, ...
            'LineStyle', '--', 'Color', 'r', 'LineWidth', 1.5);
        text(res.f_beat_theory, yl(2)*0.92, '理论值', 'Color', 'r', ...
            'FontSize', 8);
    end
    ylabel(sprintf('l=%d', res.l));
    grid on;
    if idx == num_l
        xlabel('频率 (Hz)');
    end
end
if exist('sgtitle', 'file')
    sgtitle('不同阶数花瓣模拍频频谱（红色虚线 = 理论拍频）');
end

% --- 图4: 测量误差汇总 ---
figure('Name', '测量结果汇总', 'Position', [100, 100, 800, 400]);
l_vals = [results.l];
err_vals = [results.error_pct];

bar(l_vals, err_vals, 0.6);
xlabel('花瓣模阶数 l');
ylabel('测量误差 (%)');
title(sprintf('不同阶数旋转频率测量误差 (真实值 = %.2f Hz)', f_true));
grid on;
set(gca, 'XTick', l_vals);

% --- 图5: 花瓣模旋转动画示意 (取 l=3) ---
figure('Name', '花瓣旋转示意', 'Position', [100, 100, 500, 500]);
l_show = 3;
[~, I_petal_show] = petal_mode_beam(R, TH, l_show, w0, 0);
theta_steps = 0:Omega_true*dt:2*pi;
n_steps = min(20, length(theta_steps));

for k = 1:n_steps
    theta_rot = theta_steps(k);
    % 坐标旋转法：查询点取旋转后的坐标（X,Y 为单调网格，合法）
    Xr = X*cos(theta_rot) - Y*sin(theta_rot);
    Yr = X*sin(theta_rot) + Y*cos(theta_rot);
    I_rot = interp2(X, Y, I_petal_show, Xr, Yr, 'linear', 0);
    imagesc(x_ax, y_ax, I_rot);
    axis square; axis tight;
    colormap(gca, 'hot');
    xlabel('x (mm)'); ylabel('y (mm)');
    title(sprintf('l = %d 花瓣旋转示意 (Ω = %.0f rad/s)', l_show, Omega_true));
    drawnow;
    pause(0.05);
end

%% ======================== 结果汇总表 ========================
fprintf('\n\n===== 测量结果汇总 =====\n');
fprintf('%-6s %-8s %-12s %-12s %-10s\n', ...
    'l', '花瓣数', '理论拍频(Hz)', '测量拍频(Hz)', '误差(%)');
fprintf('--------------------------------------------\n');
for idx = 1:num_l
    fprintf('%-6d %-8d %-12.2f %-12.2f %-10.3f\n', ...
        results(idx).l, results(idx).petals, ...
        results(idx).f_beat_theory, results(idx).f_peak, ...
        results(idx).error_pct);
end
fprintf('\n真实旋转频率: %.2f Hz\n', f_true);
fprintf('===== 仿真完成 =====\n');
