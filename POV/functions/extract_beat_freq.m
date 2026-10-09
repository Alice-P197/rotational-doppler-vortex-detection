function [f_axis, P1, f_peak, I_peak] = extract_beat_freq(t, I_signal, f_min, f_max)
% EXTRACT_BEAT_FREQ 从时域强度信号中提取拍频频率
%
% 对时域信号做 FFT，在指定频率范围内找到功率谱峰值
%
% 输入:
%   t       - 时间轴 (1×M)
%   I_signal - 时域强度信号 (1×M)
%   f_min   - 搜索频率下限 (Hz)
%   f_max   - 搜索频率上限 (Hz)
%
% 输出:
%   f_axis  - 频率轴单边谱 (Hz)
%   P1      - 单边幅度谱
%   f_peak  - 峰值频率 (Hz)
%   I_peak  - 峰值幅度
%
% 用法:
%   [f, P, fp, Ip] = extract_beat_freq(t, I, 1, 1000);

    % 去均值（去除直流分量）
    I_centered = I_signal - mean(I_signal);

    % FFT
    M = length(I_centered);
    dt = t(2) - t(1);
    Fs = 1 / dt;

    Y = fft(I_centered);
    P2 = abs(Y / M);             % 双边谱
    P1 = P2(1:floor(M/2)+1);      % 单边谱
    P1(2:end-1) = 2 * P1(2:end-1);
    f_axis = Fs * (0:(M/2)) / M;

    % 在指定范围内找峰值
    mask = (f_axis >= f_min) & (f_axis <= f_max);
    f_search = f_axis(mask);
    P_search = P1(mask);

    [I_peak, idx] = max(P_search);
    f_peak = f_search(idx);
end
