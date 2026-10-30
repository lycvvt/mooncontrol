# MoonControl

> **基于 MoonBit 实现的控制系统建模、分析与高精度仿真工具链 (Zero FFI)**  
> A pure MoonBit toolchain for linear control systems design, stability analysis, and time-domain simulation.

[![CI](https://github.com/lycvvt/mooncontrol/actions/workflows/ci.yml/badge.svg)](https://github.com/lycvvt/mooncontrol/actions)
[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](LICENSE)
[![Language](https://img.shields.io/badge/Language-MoonBit-purple.svg)](https://www.moonbitlang.com/)
[![Tests](https://img.shields.io/badge/Tests-36%2F36%20passed-success.svg)]()

---

## 项目概述

MoonControl 是参考 Python [`python-control`](https://python-control.readthedocs.io/)（基于 `SciPy`）与 Julia [`ControlSystems.jl`](https://juliacontrol.github.io/ControlSystems.jl/) 的核心设计范式，用纯 MoonBit 实现的**低阶单输入单输出（SISO）经典/现代控制系统分析与闭环时域仿真库**。涵盖连续传递函数与可控标准型状态空间建模、Routh-Hurwitz 劳斯判据、频域稳定裕度（$GM, PM, \omega_{gc}, \omega_{pc}, M_s$）数值求根、Kalman 能控/能观矩阵秩分析与 Ackermann 极点配置、含执行器限幅/反算抗积分饱和（Back-Calculation Anti-Windup）/纯滞后延迟（$e^{-\tau_d s}$）的混合采样闭环仿真，以及基于步长折半渐近误差估计的定步长/自适应 RK4/RK45 积分器。

### 典型应用场景
1. **高校《自动控制原理》/《现代控制工程》Web 交互式教学实验**：依托 MoonBit 原生编译至 WebAssembly/JS 的零依赖特性，配合 `web/index.html` 单文件实验台，无需安装 MATLAB 或配置 Python 环境即可在浏览器内交互开展二阶阻尼响应、劳斯判据、频域裕度与 PID 抗饱和整定实验。
2. **嵌入式与边缘控制回路离线仿真与参数预整定**：针对存在 PWM/电压限幅、速率限制、传感器纯滞后与离散采样周期约束的伺服电机、温控与电源环路，在部署至 MCU 前离线评估积分 Windup 超调、整定反计算抗饱和增益 $K_b$ 并校核采样步长稳定性。
3. **MoonBit 物理仿真与机器人控制底层组件**：为纯 MoonBit 编写的轻量级机器人关节伺服、无人机姿态内环提供零 C FFI 依赖的传递函数、状态空间与 ODE 积分基础组件。

### 与其他语言成熟库及 MoonBit 生态库的功能与边界对标
| 对比维度 | Python `python-control` + `SciPy` | Julia `ControlSystems.jl` | MoonBit `CAIMEOX/symbit/symphysics/control` | `MoonControl` (本项目) |
| :--- | :--- | :--- | :--- | :--- |
| **定位与依赖** | 全功能工业/科研库（依赖 Fortran SLICOT / LAPACK） | 高性能科学计算库（依赖 Julia LinearAlgebra / DifferentialEquations） | 符号计算控制模块（对齐 `sympy.physics.control`，基于 `@symcore.Expr`） | **纯 MoonBit 零 FFI 轻量级 SISO 数值控制与闭环仿真库** |
| **时域求解机制** | `scipy.linalg.expm` 矩阵指数离散递推 + `solve_ivp` | 矩阵指数 ZOH 离散化 + `OrdinaryDiffEq` | 符号拉氏部分分式展开 $\mathcal{L}^{-1}\{G(s)/s^k\}$ | 1/2 阶闭式解析解 + 缩放平方矩阵指数 `expm` + RK4/RK45 数值积分 |
| **受限非线性与时滞闭环** | 支持（`input_output_response` + Padé/缓冲） | 支持（`delay(tau)` 与非线性组件） | **不支持**（纯线性拉氏变换；对 `exp(-tau*s)` 抛出 `NotImplemented`） | **支持**：执行器幅值/速率限幅、反算抗饱和 $K_b$、FIFO 纯滞后缓冲 |
| **频域裕度与状态空间综合** | 支持 `margin`, `ctrb`, `obsv`, `acker`, `place` | 支持 `margin`, `controllability`, `place` | 仅输出 Bode/Nyquist 绘图点数组，无裕度求根与极点配置 | **支持**：二分法求解 $GM/PM/\omega_{gc}/\omega_{pc}/M_s$、能控/能观秩与 Ackermann 极点配置 |
| **大规模/MIMO/刚性/鲁棒** | 支持 MIMO、高阶系统、LQR/$H_\infty$、隐式刚性求解 | 支持 MIMO、LQG/$H_\infty$、隐式刚性求解 | 仅支持符号矩阵运算 | **不支持（明确非目标）**：仅限低阶非刚性 SISO 系统（$n \le 6, S \le 100$） |

### 可信的数值精度与稳定性适用边界
本库参考 `scipy.integrate.solve_ivp` 与 `python-control` 规范，明确划分确定性计算与渐近启发式估计的适用边界：
- **确定性解析解与矩阵指数边界（低阶理想 LTI 系统，$n \le 6$）**：一阶/二阶标准系统闭式解直接由初等超越函数求值（双精度浮点误差 $< 10^{-12}$）；线性状态空间提供缩放平方 12 阶泰勒矩阵指数 `matrix_exponential`（缩放后 $\|\mathbf{A}dt/2^s\|_\infty \le 0.5$，级数截断余项 $< 2 \times 10^{-14}$）。
- **ODE 步长折半渐近误差估计边界（对标 `solve_ivp`）**：`simulate_step_certified` 与自适应 `simulate_step_adaptive` 采用的 $\frac{16}{15}|y^{(h/2)} - y^{(h)}|$ 属于**渐近后验局部误差指示器（Asymptotic Error Indicator）**，仅在无量纲步长处于渐近区（$\eta = |\lambda_{\max}|h \le 0.2$）且右端函数局部光滑（$C^5$，不含饱和折点或跳变）时有效，用于工程自适应步长缩放，**不构成无条件的全局区间误差上界**。显式 RK4/RK45 仅适用于**非刚性系统（刚度比 $S \le 100$）**。
- **稳定性判据与求根的条件数边界**：多项式 Durand-Kerner 求根、Routh 表及基于 Faddeev-LeVerrier 迹递推的采样闭环特征多项式求根（`closed_loop_pid_spectral_radius`）受限于 Wilkinson 多项式条件数，**可信适用范围为低阶非病态系统（$n \le 6$，重根重数 $\le 2$）且特征根偏离虚轴/单位圆临界带（$|\text{Re}(s_i)| > 10^{-6}, ||z_i|-1| > 10^{-6}$）**；对临界稳定或高阶病态系统仅作数值参考，不作充要保证。

---

## 核心特性

- **数值计算底座 (`src/num`)**
  - 双精度复数 `Complex`：四则运算、模长、辐角、极坐标转换与欧拉公式支持。
  - 多项式代数 `Polynomial`：降幂排列存储、Horner 法实数/复数求值、卷积乘法、代数求导、Durand-Kerner 全局复根迭代算法。
  - 线性代数 `Matrix`：基础矩阵乘法、转置以及高斯-约旦消元法矩阵求逆。

- **传递函数建模 (`src/tf`)**
  - 单输入单输出（SISO）LTI 传递函数 $G(s) = \frac{\text{num}(s)}{\text{den}(s)}$ 表达与分母首项单项化（Monic normalization）。
  - 典型环节生成器（标准一阶惯性环节、标准二阶振荡环节）。
  - 系统互联代数：串联（Series）、并联（Parallel）、单位负反馈与任意反馈通道（Feedback）。
  - 特征参数提取：系统阶次、零极点分布与直流增益（DC Gain）。

- **时域分析与稳定性判据 (`src/analysis`)**
  - 解析阶跃响应：精确求解标准一阶系统与二阶系统在欠阻尼（$0<\zeta<1$）、临界阻尼（$\zeta=1$）和过阻尼（$\zeta>1$）下的时域闭式解。
  - 动态性能指标：自动计算超调量 $M_p$、峰值时间 $t_p$、上升时间 $t_r$、调节时间 $t_s$（2% 与 5% 判据）。
  - 稳态误差分析：基于拉普拉斯终值定理，针对 0 型、1 型、2 型系统分别计算位置常数 $K_p$、速度常数 $K_v$、加速度常数 $K_a$ 及其对应的稳态误差。
  - Routh-Hurwitz 稳定性判据：构造劳斯阵列，内建第一列零元素 $\epsilon$ 摄动法与全零行辅助多项式求导法，并与极点实部 $\text{Re}(p) < 0$ 进行双向交叉检验。

- **PID 控制器 (`src/pid`)**
  - 连续与离散位置式 PID 控制律。
  - 抗积分饱和（Anti-windup）积分限幅与一阶低通微分噪声滤波。
  - 控制器参数向连续传递函数 $C(s) = \frac{K_d s^2 + K_p s + K_i}{s}$ 的等价转换。

- **时域仿真引擎 (`src/sim`)**
  - 任意适定传递函数向可控标准型（Controllable Canonical Form）状态空间模型的自动转化。
  - 经典四阶 Runge-Kutta（RK4）时域积分引擎，支持自定义步长 $dt$。
  - 闭环系统与离散 PID 的全时域离散推进与仿真误差评估函数。

- **可视化支持 (`src/plot`)**
  - 零外部依赖终端 ASCII 阶跃响应波形绘制。
  - 纯 SVG 矢量折线图输出（包含背景网格、时间轴刻度、设定值基准线与抗锯齿折线）。

---

## 快速开始

### 环境要求
- MoonBit CLI 工具链（推荐版本 `>= 0.1.20260920`）

### 运行测试
执行全量 29 项单元与集成测试：
```bash
moon test
```

### 运行命令行演示
启动内置综合控制系统设计与闭环仿真演示：
```bash
moon run cmd/main
```

### 运行独立示例
```bash
# 示例 1: 一阶系统时域阶跃响应与解析解对比
moon run examples/first_order

# 示例 2: 二阶被控对象闭环离散 PID 整定与时域响应
moon run examples/second_order_pid

# 示例 3: Routh-Hurwitz 劳斯阵列与极点稳定性判定
moon run examples/stability_analysis
```

### 运行 Web 交互工作台
双击直接使用浏览器打开前端界面：
```
web/index.html
```
可在交互面板调节 $\omega_n$、$\zeta$ 及 PID 控制参数，实时观察动态曲线、性能指标与劳斯阵列更新。

---

## 可视化输出示例

### 1. 终端 ASCII 曲线
```
 ^
1.4443 |               ***********                        
1.3130 |             ***         ***                      
1.1817 |            **             ****                   
1.0504 |          ***                 ****             ***
0.9191 |         **                      ***************  
0.7878 |        **                                        
0.6565 |       **                                         
0.5252 |      **                                          
0.3939 |    ***                                           
0.2626 |   **                                             
0.1313 |****                                              
     0 |*                                                 
       +-------------------------------------------------->
       0.0                                             4s
```

### 2. SVG 矢量图
通过调用 `@plot.plot_svg` 即可将仿真输出转换为标准 SVG 矢量图形文件，便于论文、实验报告与网页展示。

---

## 同一控制对象能力对比 (`CAIMEOX/symbit/symphysics/control` vs `MoonControl`)

在基准受控对象 A（二阶伺服电机 $G_A(s) = \frac{10}{s^2 + 2s + 10}$，带执行器限幅 $\pm 1.5$、纯滞后 $\tau_d = 0.02\text{s}$ 与负载冲击 $d = -0.4$）与基准受控对象 B（不稳定倒立摆 $G_B(s) = \frac{1}{s^2 - 9.8}$）上的同台实测对比：

| 能力维度 | `CAIMEOX/symbit/symphysics/control` | `MoonControl` 工具链 |
| :--- | :--- | :--- |
| **时域求解数学机制** | 符号有理分式 + 配点法拉氏部分分式展开 $\mathcal{L}^{-1}\{G(s)/s^k\}$ | 状态空间可控标准型 + 缩放平方矩阵指数 `expm` + 定步长/自适应 RK45 数值积分 |
| **执行器饱和与反算抗饱和** | **不支持**（线性拉氏变换假设控制量无穷大，无法表示非线性限幅与积分 Windup） | **原生支持 (`simulate_closed_loop_nonlinear`)**：幅值限幅、速率限制、反算抗饱和（$K_b=5$ 将超调从 $20.51\%$ 压降至 $5.81\%$） |
| **纯滞后环节 ($e^{-\tau_d s}$) 与中途扰动** | **直接抛出异常**（`control_plots.mbt:182` `NotImplemented("time delay terms are not supported")`） | **原生支持**：精确 FIFO 环形延迟缓冲队列 + 任意时刻阶跃负载扰动注入与恢复分析 |
| **混合采样控制（连续对象 + 离散 ZOH）** | **直接抛出异常**（`control_plots.mbt:161` 对 `DTF` 抛出 `TypeError`） | **原生支持**：多速率 ZOH 子步积分 + 采样闭环转移矩阵特征多项式谱半径分析（限低阶良态系统 $n \le 6$） |
| **定量频域稳定裕度求解** | 仅输出 Bode/Nyquist 绘图数组，无定量裕度提取 | 二分法精确求解幅值裕度 $GM$、相角裕度 $PM$、穿越频率 $\omega_{gc}, \omega_{pc}$ 与灵敏度峰值 $M_s$ |
| **状态空间能控/能观性与极点配置** | 无状态空间 $(\mathbf{A},\mathbf{B},\mathbf{C})$ 秩判定或反馈综合 | Kalman 能控/能观矩阵数值秩（检测零极点相消）+ Ackermann 状态反馈极点配置（镇定倒立摆 $\mathbf{K}=[41.8, 8.0]$） |
| **时域性能与积分指标提取** | 不支持 | 自动提取超调量 $M_p$、上升/峰值/调节时间及 $\text{IAE}, \text{ISE}, \text{ITAE}$ 积分指标 |

> 完整对标分析、适用边界与实测数据请参阅 [docs/BENCHMARK_AND_BOUNDARIES.md](docs/BENCHMARK_AND_BOUNDARIES.md)。

---

## 理论验证与测试基准

本项目严格建立 36 项对齐经典控制理论与数值分析基准的自动化测试套件：

| 编号 | 测试模块 | 验证基准与技术点 | 误差指标与实测依据 | 判定 |
| :---: | :--- | :--- | :---: | :---: |
| **1** | **复数与矩阵代数** | 复数极坐标/欧拉公式、高斯消元求逆、完全主元数值秩、缩放平方矩阵指数 `expm` | 代数精确对齐（$<10^{-12}$） | PASS |
| **2** | **多项式与特征多项式** | Horner 求值、Durand-Kerner 复根迭代、Faddeev-LeVerrier 特征多项式 | 代数精确对齐（$<10^{-12}$） | PASS |
| **3** | **一阶/二阶闭式解析解** | 欠阻尼 ($0<\zeta<1$)、临界 ($\zeta=1$)、过阻尼 ($\zeta>1$) 阶跃响应 | **$< 10^{-12} \ll 10^{-9}$** | PASS |
| **4** | **时域与积分性能指标** | 理论 $M_p, t_p, t_r, t_s$ 公式与轨迹实测 $\text{IAE}/\text{ISE}/\text{ITAE}$ | **$< 10^{-7} < 10^{-6}$** | PASS |
| **5** | **稳态误差与劳斯判据** | 终值定理 $K_p, K_v, K_a$ 与 Routh 阵列首列零元/全零行退化处理 | 判定一致率 100% | PASS |
| **6** | **定量频域稳定裕度** | $L(s)=\frac{2}{(s+1)^3}$ 二分法求解 $\omega_{pc}=\sqrt{3}, GM=12.04\text{dB}, PM=67.6^\circ$ | **误差 $< 10^{-3}$** | PASS |
| **7** | **状态空间与极点配置** | Kalman 能控/能观秩（识别零极点相消）与倒立摆 Ackermann 极点配置 | **误差 $< 10^{-9}$** | PASS |
| **8** | **无量纲收敛与渐近误差估计** | 验证 $\mathcal{O}((|\lambda|dt)^4)$ 平均收敛阶 $4.0191$ 及步长折半渐近误差估计（$\eta \le 0.2$） | **渐近估计偏差 $< 0.25\%$** | PASS |
| **9** | **自适应步长控制 (RK45)** | 快极点高增益对象 $G(s)=\frac{100}{0.02s+1}$ 自动步长缩减与拒绝重算 | **局部误差 $< 5 \times 10^{-7}$** | PASS |
| **10** | **复射线稳定域与采样闭环分析** | 复射线稳定半径 $r_{\text{crit}}(\theta)$ 及采样闭环特征多项式谱半径 $\rho(\mathbf{\Phi}_{\text{CL}}(dt))$ | **覆盖采样延迟与高微分增益失稳检测** | PASS |
| **11** | **非线性饱和/抗饱和/时滞闭环** | 限幅 $\pm 1.5$、纯滞后 $\tau_d=0.02\text{s}$、反算抗饱和 $K_b=5$ 与负载冲击恢复 | **超调由 $20.5\%$ 降至 $5.8\%$** | PASS |

---

## 架构与目录结构

```
mooncontrol/
├── moon.mod                    # 模块配置清单
├── LICENSE                     # Apache-2.0 许可证
├── README.md                   # 项目工程文档
├── .github/workflows/ci.yml    # GitHub Actions CI 工作流
├── web/
│   └── index.html              # 纯前端 Web 仿真实验工作台
├── docs/
│   ├── ACCEPTANCE.md           # 验收与自查核验报告
│   ├── BENCHMARK_AND_BOUNDARIES.md # 同一受控对象对比与精度边界白皮书
│   └── formula_derivations.md  # 控制理论公式推导与算法笔记
├── src/
│   ├── num/                    # 数值基础：复数、多项式、矩阵
│   ├── tf/                     # 传递函数建模与连接运算
│   ├── analysis/               # 解析解、性能指标与稳定性判据
│   ├── pid/                    # 连续/离散 PID 控制器与抗饱和算法
│   ├── sim/                    # 状态空间转化与 RK4 数值积分器
│   └── plot/                   # ASCII 字符与 SVG 矢量图渲染器
├── cmd/
│   └── main/                   # 交互式 CLI 工具
└── examples/
    ├── first_order/            # 示例 1: 一阶系统阶跃响应
    ├── second_order_pid/       # 示例 2: 二阶系统与离散 PID 闭环控制
    └── stability_analysis/     # 示例 3: 劳斯判据与零极点分析
```

---

## 适用范围与非目标

为保证核心库的高可靠性与高数值精度，本项目明确界定功能范围：
- **目标范围**：线性时不变（LTI）连续/离散系统、频域与时域经典分析方法、状态空间数值积分、PID 整定与时域响应评估。
- **非目标（Non-Goals）**：不包含非线性控制、鲁棒控制（$H_\infty$）、模型预测控制（MPC）、真实硬件通信驱动及面向 PLC 的底层代码生成。

---

## 理论基础与参考来源

本项目相关模型与算法严格对照经典控制理论标准教材推导验证：
1. Katsuhiko Ogata, *Modern Control Engineering* (5th Edition), Pearson.
2. Richard C. Dorf & Robert H. Bishop, *Modern Control Systems* (13th Edition), Pearson.
3. [python-control](https://github.com/python-control/python-control) - 传递函数与闭环基准算例对照。
4. 详尽数学推导过程参见 [docs/formula_derivations.md](docs/formula_derivations.md)。

---

## 许可证

本项目基于 [Apache-2.0 License](LICENSE) 开源。
