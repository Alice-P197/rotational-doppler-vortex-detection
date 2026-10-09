%% 不同拓扑荷数 l 在相同入射条件下的旋转测速对比
%  模型与参数完全基于《涡旋光任意入射条件下的旋转物体转速探测》
%  （红外与激光工程 2021, 50(9): 20210451）复现脚本 reproduce_main.m
%
%  相同条件（所有 l 一致）:
%    入射几何 γ、φ、d、r，旋转频率 f_rot，仿真时长 T、采样率 Fs
%  对比维度:
%    1) 任意入射频移分布 Δf(θ)（公式10）随 l 的整体缩放
%    2) 特征点频移 f1~f4（θ=0,π/2,π,3π/2）与公式(12) f_mod 测速精度
%    3) 频移带宽随 l 变化（∝l）与相对带宽（与 l 无关）
%    4) 频率分辨率限制: 特征频移间距 vs δf=1/T —— 判定最小可分辨 l
%    5) 散射点频谱仿真展宽对比
%
%  输出: results_l_compare/ 下 4 张图 + 控制台汇总表

clear; clc; close all;
scriptDir = fileparts(mfilename('fullpath'));

outDir = fullfile(scriptDir, 'results_l_compare');
if ~exist(outDir, 'dir'), mkdir(outDir); end

%% ============ 相同条件（文献图3参数） ============
f_rot = 50;              % 物体旋转频率 (Hz)
Omega = 2*pi*f_rot;
r     = 4e-3;            % 涡旋光半径 (m)
d     = 1e-3;            % 光斑中心横向位移 (m)
gamma = 30*pi/180;       % 倾斜角
phi   = 30*pi/180;       % 光斑-旋转中心连线夹角

l_list = [1 2 4 6 10 18 30];   % 对比的拓扑荷数

% 频谱仿真参数（文献 [4] 散射点模型）
T_sim = 2.0;             % 采样时长 s
Fs    = 50e3;            % 采样率 Hz
Np    = 400;             % 散射点数
df_res = 1/T_sim;        % FFT 频率分辨率 δf

theta   = linspace(0, 2*pi, 2001);
th_feat = [0 pi/2 pi 3*pi/2];

fprintf('===== 相同入射条件下不同拓扑荷数测速对比 =====\n');
fprintf('相同条件: γ=30°, φ=30°, d=1mm, r=4mm, f=%.0fHz, T=%.1fs, δf=%.1fHz\n\n', ...
    f_rot, T_sim, df_res);

%% ============ 逐 l 计算 ============
num_l = numel(l_list);
results = struct();

for idx = 1:num_l
    l = l_list(idx);
    base = l * Omega / pi;   % 叠加态对准频移 lΩ/π

    % 公式(10) 任意入射频移分布
    Df = base .* ( sin(theta).^2 / cos(gamma) + cos(gamma) .* cos(theta).^2 ...
                 + (d/r) .* (sin(theta).*cos(phi) - cos(theta).*sin(phi).*cos(gamma)) );

    % 特征点频移 f1~f4
    f1 = base * (cos(gamma) - (d/r)*sin(phi)*cos(gamma));   % θ=0
    f2 = base * (1/cos(gamma) + (d/r)*cos(phi));            % θ=π/2
    f3 = base * (cos(gamma) + (d/r)*sin(phi)*cos(gamma));   % θ=π
    f4 = base * (1/cos(gamma) - (d/r)*cos(phi));            % θ=3π/2

    % 公式(12) 几何平均提取 f_mod → 转速
    f_mod  = sqrt( ((f1+f3)/2) * ((f2+f4)/2) );
    f_est  = f_mod / (2*l);
    err_pct = abs(f_est - f_rot) / f_rot * 100;

    % 带宽
    bw     = max(Df) - min(Df);
    rel_bw = bw / base;

    % 可分辨性: 特征频移最小间距 vs 频率分辨率
    feat_sorted = sort([f1 f2 f3 f4]);
    gap_min     = min(diff(feat_sorted));
    resolvable  = gap_min > df_res;

    % 散射点频谱仿真（用于展示展宽）
    th_k = (0:Np-1)/Np * 2*pi;
    Df_k = base .* ( sin(th_k).^2/cos(gamma) + cos(gamma).*cos(th_k).^2 ...
                   + (d/r).*(sin(th_k).*cos(phi) - cos(th_k).*sin(phi).*cos(gamma)) );
    tvec = (0:round(T_sim*Fs)-1) / Fs;
    s_sim = sum(exp(1i*2*pi*Df_k(:)*tvec), 1);
    f_ax  = linspace(0, Fs/2, numel(tvec)/2+1);
    P_sim = abs(fft(s_sim(:))) / numel(s_sim);
    P_sim = P_sim(1:numel(f_ax));

    results(idx).l        = l;
    results(idx).base     = base;
    results(idx).Df       = Df;
    results(idx).f_feat   = [f1 f2 f3 f4];
    results(idx).f_mod    = f_mod;
    results(idx).f_est    = f_est;
    results(idx).err_pct  = err_pct;
    results(idx).bw       = bw;
    results(idx).rel_bw   = rel_bw;
    results(idx).gap_min  = gap_min;
    results(idx).resolvable = resolvable;
    results(idx).f_ax     = f_ax;
    results(idx).P_sim    = P_sim;

    if resolvable, tag = '可分辨'; else, tag = '不可分辨'; end
    fprintf('l=%2d: 基准 %6.1f Hz | f1..f4 = %6.1f/%6.1f/%6.1f/%6.1f Hz | f_mod=%6.1f Hz | 转速误差 %6.3f%% | 带宽 %6.1f Hz | 特征间距 %4.1f Hz (>δf=%0.1f: %s)\n', ...
        l, base, f1, f2, f3, f4, f_mod, err_pct, bw, gap_min, df_res, tag);
end

%% ============ 图1: 频移分布随 l 缩放 ============
fig = figure('Position', [100 100 780 480], 'Color', 'w');
cols = lines(4);
plot_l = [1 4 10 30];
lg = gobjects(0);
for k = 1:numel(plot_l)
    idx = find([results.l] == plot_l(k));
    plot(theta, results(idx).Df, '-', 'Color', cols(k,:), 'LineWidth', 1.7); hold on;
    lg(k) = plot(NaN, NaN, '-', 'Color', cols(k,:), 'LineWidth', 1.7);
end
xlim([0 2*pi]); xticks(0:pi/2:2*pi);
xticklabels({'0','π/2','π','3π/2','2π'});
grid on; box on;
xlabel('\theta / rad'); ylabel('Rotational frequency / Hz');
title(sprintf('任意入射频移分布 \\Deltaf(\\theta) 随 l 缩放 (\\gamma=30°, \\phi=30°, d=1mm, f=%.0fHz)', f_rot));
legend(lg, arrayfun(@(x) sprintf('l=%d', x), plot_l, 'UniformOutput', false), 'Location', 'northwest');
saveas(fig, fullfile(outDir, 'fig_delf_vs_l.png')); close(fig);

%% ============ 图2: 特征频移、带宽、相对带宽 ============
fig = figure('Position', [100 100 780 480], 'Color', 'w');
lv = [results.l];
feat_all = reshape([results.f_feat], 4, num_l)';
plot(lv, feat_all, '-o', 'LineWidth', 1.5, 'MarkerSize', 6); hold on;
plot(lv, [results.bw], 'k--', 'LineWidth', 1.8);
xlabel('拓扑荷数 l'); ylabel('频率 (Hz)');
legend({'f_1 (θ=0)', 'f_2 (θ=π/2)', 'f_3 (θ=π)', 'f_4 (θ=3π/2)', '频移带宽 Δf_{max}-Δf_{min}'}, ...
    'Location', 'northwest');
title('特征点频移与带宽随 l 线性增长（相对带宽与 l 无关）');
grid on; box on;
set(gca, 'XTick', lv);
saveas(fig, fullfile(outDir, 'fig_feature_vs_l.png')); close(fig);

%% ============ 图3: 频谱展宽对比 (l=1, 4, 18) ============
spec_l = [1 4 18];
fig = figure('Position', [100 100 1100 700], 'Color', 'w');
for k = 1:numel(spec_l)
    idx = find([results.l] == spec_l(k));
    res = results(idx);
    mask = (res.f_ax >= res.base*0.3) & (res.f_ax <= res.base*1.7);
    subplot(2, 2, k);
    plot(res.f_ax(mask), 20*log10(res.P_sim(mask)/max(res.P_sim)+eps), 'b-', 'LineWidth', 1.3); hold on;
    yl = ylim;
    xline(min(res.Df), 'r--', 'LineWidth', 1.2);
    xline(max(res.Df), 'r--', 'LineWidth', 1.2);
    xlim([res.base*0.3 res.base*1.7]); ylim([-60 2]);
    grid on; box on;
    xlabel('Frequency / Hz'); ylabel('Relative amplitude / dB');
    title(sprintf('l=%d: 展宽 [%.0f, %.0f] Hz', spec_l(k), min(res.Df), max(res.Df)));
end
if exist('sgtitle', 'file')
    sgtitle('散射点频谱展宽对比（红虚线=理论频移边界）');
end
saveas(fig, fullfile(outDir, 'fig_spectra_vs_l.png')); close(fig);

%% ============ 图4: 可分辨性 vs l ============
fig = figure('Position', [100 100 780 480], 'Color', 'w');
yyaxis left;
semilogy(lv, [results.gap_min], 'bo-', 'LineWidth', 1.6, 'MarkerSize', 7, 'MarkerFaceColor', 'b');
yline(df_res, 'r--', 'LineWidth', 1.6);
ylabel('特征频移最小间距 (Hz, 对数)');
yyaxis right;
plot(lv, [results.rel_bw], 'gs-', 'LineWidth', 1.5, 'MarkerSize', 6);
ylabel('相对带宽 (Δf_{max}-Δf_{min})/lΩ/π');
xlabel('拓扑荷数 l');
title(sprintf('可分辨性判据: 特征间距 > 频率分辨率 δf=%.1f Hz (T=%.1fs)', df_res, T_sim));
legend({'最小特征间距', sprintf('δf = 1/T = %.1f Hz', df_res), '相对带宽'}, 'Location', 'southwest');
grid on; box on;
set(gca, 'XTick', lv);
saveas(fig, fullfile(outDir, 'fig_resolvable_vs_l.png')); close(fig);

%% ============ 汇总表 ============
fprintf('\n===== 汇总表 =====\n');
fprintf('%-4s %-9s %-20s %-10s %-10s %-14s %-14s\n', ...
    'l', '基准lΩ/π', '特征频移 f1/f2/f3/f4 (Hz)', 'f_mod', '误差%', '带宽(Hz)', '可分辨');
fprintf('--------------------------------------------------------------------------------\n');
for idx = 1:num_l
    fv = results(idx).f_feat;
    if results(idx).resolvable, tag2 = '是'; else, tag2 = '否'; end
    fprintf('%-4d %-9.0f %-20s %-10.1f %-10.3f %-14.1f %-14s\n', ...
        results(idx).l, results(idx).base, ...
        sprintf('%.0f/%.0f/%.0f/%.0f', fv(1), fv(2), fv(3), fv(4)), ...
        results(idx).f_mod, results(idx).err_pct, results(idx).bw, tag2);
end
fprintf('\n图片已保存至 %s\n', outDir);
fprintf('===== 对比完成 =====\n');
