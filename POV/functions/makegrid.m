function [X, Y, R, TH] = makegrid(N, L)
% MAKEGRID 生成二维直角坐标与极坐标网格
%
% 输入:
%   N - 采样点数（正方形，N×N）
%   L - 视场大小（m，正方形边长）
%
% 输出:
%   X, Y  - 直角坐标网格 (N×N)
%   R     - 径向坐标 (N×N)
%   TH    - 角向坐标 (N×N)，范围 [-π, π]
%
% 用法:
%   [X, Y, R, TH] = makegrid(512, 10e-3);

    dx = L / N;
    x = linspace(-L/2 + dx/2, L/2 - dx/2, N);
    [X, Y] = meshgrid(x, x);
    R = sqrt(X.^2 + Y.^2);
    TH = angle(X + 1i*Y);
end
