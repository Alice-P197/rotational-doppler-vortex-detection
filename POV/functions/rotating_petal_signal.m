function [t, I_signal] = rotating_petal_signal(R, TH, l, w0, Omega, T, dt, seed)
% ROTATING_PETAL_SIGNAL 模拟旋转花瓣模光束照射旋转粗糙表面的时域强度信号
%
% 物理模型:
%   花瓣模强度 I(r,θ) ∝ cos²(lθ)，照射旋转物体（表面反射率 σ）
%   探测器信号: I_det(t) = ∬ I(r, θ-Ωt)·σ(r,θ) r dr dθ
%   表面含方位标记(cos(2lθ)调制) → 产生旋转多普勒拍频
%       f_beat = 2l·Ω / (2π) = l·Ω / π
%
% 实现: 角向傅里叶法（快速、数值稳定）
%   I_det(t) = Σ_n c_n·exp(i·n·Ω·t)
%   花瓣模强度只有 0 与 ±2l 阶角向谐波非零，
%   因此只需算 3 个傅里叶系数即可解析得到拍频，
%   避免了逐时刻插值整张图的巨大开销。
%
% 输入:
%   R, TH  - 极坐标网格 (N×N)，由 makegrid 生成
%   l      - 花瓣模阶数（花瓣数 = 2l）
%   w0     - 束腰半径 (m)
%   Omega  - 物体角速度 (rad/s)
%   T      - 总仿真时间 (s)
%   dt     - 时间步长 (s)
%   seed   - 随机数种子（默认 1）
%
% 输出:
%   t       - 时间轴 (1×M)
%   I_signal - 时域强度信号 (1×M)
%
% 用法:
%   [t, I] = rotating_petal_signal(R, TH, 3, 1e-3, 2*pi*50, 0.2, 1e-5, 42);

    if nargin < 8
        seed = 1;
    end
    rng(seed);  % 固定随机种子保证可复现

    % 1. 花瓣模光束强度分布
    [~, I_petal] = petal_mode_beam(R, TH, l, w0, 0);

    % 2. 旋转物体表面反射率
    %    粗糙随机场 + 确定性方位标记（保证 ±2l 阶拍频清晰稳定）
    N = size(R, 1);
    phi0 = 0.3;
    sigma = 0.3 + 0.08*cos(2*l*TH + phi0) + 0.1*randn(N);
    sigma(sigma < 0) = 0;

    % 3. 角向离散化（一维），径向用 r 加权累加
    K = 2048;                     % 角向采样数（> 4l 即可，2048 富余）
    dtheta = 2*pi / K;
    idx_ang = mod(round((TH + pi) / dtheta), K);   % 0..K-1

    w = abs(R);                    % 面积权重 ∝ r·dr
    I_ang = accumarray(idx_ang(:)+1, I_petal(:).*w(:), [K 1]);  % K×1
    S_ang = accumarray(idx_ang(:)+1, sigma(:).*w(:),   [K 1]);  % K×1

    % 4. FFT 角向频谱
    I_F = fft(I_ang);              % I_n
    S_F = fft(S_ang);              % σ_n

    % 5. 时间轴
    t = 0:dt:T;
    M = length(t);

    % 6. 解析拍频信号
    %    I_det(t) = c0 + 2·Re[ c_{2l}·exp(-i·2l·Ω·t) ]
    %    系数: c_n = (2π/K)·I_n·conj(σ_n)
    c0  = (2*pi/K) * I_F(1) * conj(S_F(1));
    c2l = (2*pi/K) * I_F(2*l+1) * conj(S_F(2*l+1));

    I_signal = real(c0) + 2*real(c2l * exp(-1i*2*l*Omega*t));

    % 7. 归一化（去直流、归一化到均值 1 附近）
    I_signal = I_signal / mean(I_signal);
end
