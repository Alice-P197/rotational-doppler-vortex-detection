function [E, I] = petal_mode_beam(R, TH, l, w0, p)
% PETAL_MODE_BEAM 生成花瓣模光束（Petal-mode beam）复振幅
%
% 花瓣模 = 两个相反 OAM 阶数的 LG 光束复振幅叠加（对齐笔记/GitHub 定义）
%   E(r,θ) = LG_{p,+l}(r,θ) + LG_{p,-l}(r,θ)
%
% LG 光束（源平面 z=0，束腰 w0）:
%   LG_{p,l} = (√2·r/w0)^|l| · exp(-r²/w0²) · L_p^|l|(2r²/w0²)
%              · exp(-i·l·θ)     （z=0 时球面相位=1、古伊相位=0）
%
% 叠加后: E ∝ (√2·r/w0)^|l| · exp(-r²/w0²) · L_p^|l|(2r²/w0²) · cos(l·θ)
% 强度呈 2l 个花瓣，中心为暗斑（相位奇点，(r)^|l| → 0）
%
% 输入:
%   R, TH - 极坐标网格 (N×N)，由 makegrid 生成
%   l     - OAM 阶数（拓扑荷数），花瓣数 = 2l
%   w0    - 束腰半径 (m)
%   p     - 径向阶数（默认 0）
%
% 输出:
%   E - 复振幅 (N×N)
%   I - 强度 |E|² (N×N)
%
% 用法:
%   [E, I] = petal_mode_beam(R, TH, 3, 1e-3);

    if nargin < 5
        p = 0;
    end

    % LG 径向振幅（z=0 处 w_z = w0）
    % (√2·r/w0)^|l| · exp(-r²/w0²) · L_p^|l|(2r²/w0²)
    radial = (sqrt(2)*R/w0).^abs(l) .* exp(-R.^2/w0^2) ...
             .* laguerre(p, abs(l), 2*(R/w0).^2);

    % 角向叠加: exp(-i·l·θ) + exp(+i·l·θ) = 2·cos(l·θ)
    E = radial .* (exp(-1i*l*TH) + exp(1i*l*TH));

    % 归一化（与笔记 LG_beam.m 一致）
    E = E / sqrt(sum(abs(E(:)).^2));

    if nargout > 1
        I = abs(E).^2;
    end
end

function result = laguerre(p, l, x)
% 广义拉盖尔多项式 L_p^l(x)（与笔记 LG_beam.m 完全一致的递归定义）
    if p == 0
        result = ones(size(x));
    elseif p == 1
        result = 1 + abs(l) - x;
    else
        result = (1/p) * ((2*p + l - 1 - x).*laguerre(p-1, abs(l), x) ...
                 - (p + l - 1)*laguerre(p-2, abs(l), x));
    end
end
