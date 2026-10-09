# 涡旋光旋转多普勒测速（Rotational Doppler Velocimetry with Vortex Beams）

基于**光学旋转多普勒效应（RDE）**的旋转物体转速探测 MATLAB 仿真，包含两部分：

1. **`POV/` — 花瓣模（Petal-mode）旋转多普勒仿真**
   用叠加态拉盖尔-高斯光束（花瓣模）照射旋转平板，仿真反射光强的拍频信号，从频谱提取旋转频率。包含花瓣模光束生成、角向傅里叶法快速拍频仿真、拍频提取三个核心模块。

2. **`reproduce_ding2021/` — 文献复现**
   复现丁友等《涡旋光任意入射条件下的旋转物体转速探测》（红外与激光工程 2021, 50(9): 20210451）的理论频移分布（公式 10）、特征点提取方法（公式 11/12）与频谱展宽现象。

---

## 物理原理

涡旋光携带轨道角动量，波前为螺旋面，坡印廷矢量相对光轴倾斜（倾斜角 α 满足 $\sin\alpha = l/(kr)$）。照射旋转平板时，表面散射点的旋转线速度在坡印廷矢量方向产生非零投影，形成旋转多普勒频移。

**对准入射**（光轴与转轴重合），叠加态 ±l 光束的频移为：

$$\Delta f = \frac{l\Omega}{\pi}$$

**任意入射**（倾斜角 γ、横向位移 d、夹角 φ），频移随散射点方位角 θ 分布（文献公式 10）：

$$\Delta f(\theta)=\frac{l\Omega}{\pi}\left[\frac{\sin^{2}\theta}{\cos\gamma}+\cos\gamma\cos^{2}\theta+\frac{d}{r}\big(\sin\theta\cos\varphi-\cos\theta\sin\varphi\cos\gamma\big)\right]$$

**频率提取**（文献公式 12，几何平均形式）：选取 θ = 0、π/2、π、3π/2 四个特征点，对偶点平均后 d、φ 项精确抵消：

$$f_{mod}=\sqrt{\frac{\Delta f_0+\Delta f_\pi}{2}\cdot\frac{\Delta f_{\pi/2}+\Delta f_{3\pi/2}}{2}}=\frac{l\Omega}{\pi}$$

物体旋转频率 $f = f_{mod}/(2l)$。**注意**：文献原文公式(12)排版丢掉了根号，本仓库采用几何平均形式（√(2012×1801)=1903.6 Hz，与文献实验值 1903.58 Hz 吻合验证）。

---

## 仓库结构

```
rotational-doppler-vortex-detection/
├── POV/                          # 花瓣模旋转多普勒仿真
│   ├── main_rotational_doppler.m # 主脚本：5 张图（花瓣分布/时域/频谱/误差/动画）
│   └── functions/
│       ├── makegrid.m            # 极坐标网格生成
│       ├── petal_mode_beam.m     # 花瓣模复振幅（LG(+l)+LG(-l)，含 (√2r/w)^|l| 暗心因子）
│       ├── rotating_petal_signal.m # 角向傅里叶法快速拍频信号生成
│       └── extract_beat_freq.m   # FFT 拍频提取
└── reproduce_ding2021/           # 丁友 2021 文献复现
    ├── reproduce_main.m          # 主脚本：公式10 曲线/特征点/频谱展宽/四工况对比
    ├── create_notion_page.py     # Notion API 建页脚本（可选）
    └── results/                  # 复现结果图（PNG）
```

---

## 使用方法

### 花瓣模旋转多普勒仿真（POV）

```matlab
cd POV
main_rotational_doppler
```

输出：l = 1~5 花瓣模的拍频测量，理论拍频 $f_{beat} = 2l\cdot f$，反推转速误差 < 0.01%。

### 文献复现（reproduce_ding2021）

```matlab
cd reproduce_ding2021
reproduce_main
```

输出 5 张图到 `results/`，并在命令行打印特征点频移、提取转速与误差。

---

## 复现结果

| 项目 | 结果 |
|---|---|
| 基准频移 lΩ/π（γ=30°, φ=30°, d=1mm, r=4mm, l=±18, f=50Hz） | 1800.0000 Hz |
| 特征点提取 f_mod | 1800.0000 Hz |
| 反推转速 | 50.0000 Hz（真值 50 Hz，**误差 0.0000%**） |
| 任意入射频移范围 | [1302.71, 2481.44] Hz |
| 频谱展宽区间 | [1303, 2481] Hz（与理论边界一致） |
| 花瓣模中心强度（l=1 / l=3） | 0.0008 / 0.0000（暗斑，相位奇点） |

复现结果图（`reproduce_ding2021/results/`）：

| 图 | 内容 |
|---|---|
| fig3_theory_curve.png | 任意入射频移分布 Δf(θ)（文献图3 复现） |
| fig_special_points.png | 四个特征点与基准频移 |
| fig_spectrum_aligned.png | 对准入射单峰频谱 |
| fig_spectrum_general.png | 任意入射频谱展宽 |
| fig_cases_compare.png | 四种入射工况频移对比 |

---

## 环境

- MATLAB R2025（R2018b+ 均可，脚本对 xline/sgtitle 做了旧版兼容）
- 无需额外工具箱，全部使用基础函数

## 引用

复现文献：

> 丁友, 丁源圣, 邱松, 刘通, 任元. 涡旋光任意入射条件下的旋转物体转速探测(特邀)[J]. 红外与激光工程, 2021, 50(9): 20210451. DOI: 10.3788/IRLA20210451

相关背景：

> Allen L, et al. Orbital angular momentum of light and the transformation of Laguerre-Gaussian laser modes [J]. Physical Review A, 1992, 45: 8185-8189.
>
> Lavery M P J, et al. Detection of a spinning object using light's orbital angular momentum [J]. Science, 2013, 341(6145): 537-540.

## 许可

仅用于学术研究与学习交流，引用请注明出处。
