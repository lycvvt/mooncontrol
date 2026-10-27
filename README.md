# MoonControl

> **基于 MoonBit 实现的控制系统建模、分析与高精度仿真工具链 (Zero FFI)**  
> A pure MoonBit toolchain for linear control systems design, stability analysis, and time-domain simulation.

[![CI](https://github.com/lycvvt/mooncontrol/actions/workflows/ci.yml/badge.svg)](https://github.com/lycvvt/mooncontrol/actions)
[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](LICENSE)
[![Language](https://img.shields.io/badge/Language-MoonBit-purple.svg)](https://www.moonbitlang.com/)
[![Tests](https://img.shields.io/badge/Tests-35%2F35%20passed-success.svg)]()

---

## 项目概述

MoonControl 是面向 MoonBit 生态构建的**状态空间常微分方程（ODE）数值闭环仿真与定量控制综合工具链**，涵盖连续 LTI 传递函数与状态空间建模、Kalman 能控/能观性秩判定、Ackermann 状态反馈极点配置、定量频域稳定裕度（$GM, PM, \omega_{gc}, \omega_{pc}, M_s$）求解、支持执行器限幅/反算抗饱和（Back-Calculation Anti-Windup）/纯滞后延迟（$e^{-\tau_d s}$）的混合采样闭环仿真，以及带 Richardson 后验误差证书与自适应步长控制（Adaptive RK45）的数值积分引擎。

### 生态定位与同 `CAIMEOX/symbit/symphysics/control` 的核心差异
在 Mooncakes 官方中心中，已有符号物理控制库 [`CAIMEOX/symbit/symphysics/control`](https://mooncakes.io/docs/CAIMEOX/symbit/symphysics/control)（对齐 Python `sympy.physics.control`）实现了基于 `@symcore.Expr` 的符号传递函数代数、符号 Routh 表以及基于**拉氏逆变换部分分式展开**（`control_partial_fraction_terms`）的线性开环响应采样。

MoonControl 在重合的线性传递函数基础层之上，针对**符号拉氏变换在原理上无法处理的工程闭环仿真与定量控制综合问题**建立了核心差异：
- **非线性与受限执行器闭环时域仿真**：拉氏部分分式展开仅对严格线性无约束系统成立；MoonControl 原生仿真执行器幅值饱和 $[u_{\min}, u_{\max}]$、速率限制 $|\dot{u}| \le r_{\max}$ 与反算抗积分饱和（Back-Calculation $K_b$）。
- **纯滞后时间延迟 ($e^{-\tau_d s}$) 与时域外部扰动注入 $d(t)$**：`symbit` 对含 `exp(-tau*s)` 的系统直接抛出 `NotImplemented("time delay terms are not supported")`；MoonControl 通过精确环形延迟队列与子步积分原生支持时滞闭环与中途负载冲击测试。
- **混合采样控制（连续对象 + 离散 ZOH 控制器）与采样闭环谱半径判据**：`symbit` 对离散传递函数绘图直接抛出 `TypeError`；MoonControl 支持零阶保持（ZOH）多速率子步积分，并通过增广采样闭环转移矩阵谱半径 $\rho(\mathbf{\Phi}_{\text{CL}}(T_s)) < 1$ 严格判定数字闭环稳定性。
- **定量频域稳定裕度与状态空间极点配置**：`symbit` 仅输出绘图采样点数组；MoonControl 二分法精确求解 $GM, PM, \omega_{gc}, \omega_{pc}, M_s$，并提供状态空间能控/能观矩阵秩检测与 Ackermann 状态反馈极点配置。

### 数值精度保证与稳定性数学边界（无反例严格界定）
- **无量纲渐近截断误差界与 Richardson 后验误差证书**：
  - 废除与时间尺度/增益脱节的固定秒数区间承诺，严格按**无量纲步长** $\eta = |\lambda_{\max}| dt$ 与模态留数 $\sum |c_i|$ 度量全局误差 $\mathcal{O}(\sum |c_i| \eta^4)$；
  - 内置步长折半 **Richardson 外推后验误差证书** $\hat{E}_{\max}(dt) = \max_k \frac{16}{15} |y^{(dt/2)}(t_k) - y^{(dt)}(t_k)|$（在基准伺服系统上与真实解析误差相对偏差 $< 0.25\%$）及**自适应步长求解器 (`simulate_step_adaptive`)**，动态缩放步长 $h_k$ 以严格满足 $(\text{atol}, \text{rtol})$ 容差。
- **开环二维复射线稳定域与闭环采样转移矩阵谱半径充要判据**：
  - *开环积分稳定界*：仅当全体极点满足严格 Hurwitz 条件 $\text{Re}(\lambda_i) < 0$ 时存在正稳定步长（开环不稳定系统强制返回 `max_stable_dt = 0.0`）；沿各极点真实复射线方向二分求解 $|R(r e^{j\theta_i})| = 1$ 得到精确 $r_{\text{crit}}(\theta_i)$，取 $dt_{\max,\text{open}} = \min_i \frac{r_{\text{crit}}(\arg \lambda_i)}{|\lambda_i|}$；
  - *采样 PID 闭环充要稳定界*：构造 $(n+2)$ 维增广状态 $\mathbf{z}_k = [\mathbf{x}_k^T, I_{k-1}, e_{k-1}]^T$ 的采样转移矩阵 $\mathbf{\Phi}_{\text{CL}}(dt)$，利用 Gelfand 幂次定理计算谱半径 $\rho(\mathbf{\Phi}_{\text{CL}}(dt)) < 1.0$，彻底消除大采样周期 ZOH 失稳与未滤波高微分增益（$c b K_d > 1$）代数环失稳反例。

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
| **时域求解数学机制** | 符号有理分式 + 配点法拉氏部分分式展开 $\mathcal{L}^{-1}\{G(s)/s^k\}$ | 状态空间可控标准型 + 定步长/自适应 RK45 微分方程数值积分 + Richardson 后验误差认证 |
| **执行器饱和与反算抗饱和** | **不支持**（线性拉氏变换假设控制量无穷大，无法表示非线性限幅与积分 Windup） | **原生支持 (`simulate_closed_loop_nonlinear`)**：幅值限幅、速率限制、反算抗饱和（$K_b=5$ 将超调从 $20.51\%$ 压降至 $5.81\%$） |
| **纯滞后环节 ($e^{-\tau_d s}$) 与中途扰动** | **直接抛出异常**（`control_plots.mbt:182` `NotImplemented("time delay terms are not supported")`） | **原生支持**：精确 FIFO 环形延迟缓冲队列 + 任意时刻阶跃负载扰动注入与恢复分析 |
| **混合采样控制（连续对象 + 离散 ZOH）** | **直接抛出异常**（`control_plots.mbt:161` 对 `DTF` 抛出 `TypeError`） | **原生支持**：多速率 ZOH 子步积分 + 采样闭环转移矩阵谱半径充要判据 $\rho(\mathbf{\Phi}_{\text{CL}}(T_s)) < 1$ |
| **定量频域稳定裕度求解** | 仅输出 Bode/Nyquist 绘图数组，无定量裕度提取 | 二分法精确求解幅值裕度 $GM$、相角裕度 $PM$、穿越频率 $\omega_{gc}, \omega_{pc}$ 与灵敏度峰值 $M_s$ |
| **状态空间能控/能观性与极点配置** | 无状态空间 $(\mathbf{A},\mathbf{B},\mathbf{C})$ 秩判定或反馈综合 | Kalman 能控/能观矩阵数值秩（检测零极点相消）+ Ackermann 状态反馈极点配置（镇定倒立摆 $\mathbf{K}=[41.8, 8.0]$） |
| **时域性能与积分指标提取** | 不支持 | 自动提取超调量 $M_p$、上升/峰值/调节时间及 $\text{IAE}, \text{ISE}, \text{ITAE}$ 积分指标 |

> 完整源码级机制对比、数学反例剖析与实测数据请参阅 [docs/BENCHMARK_AND_BOUNDARIES.md](docs/BENCHMARK_AND_BOUNDARIES.md)。

---

## 理论验证与测试基准

本项目严格建立 35 项对齐经典控制理论与数值分析基准的自动化测试套件：

| 编号 | 测试模块 | 验证基准与技术点 | 误差指标与实测依据 | 判定 |
| :---: | :--- | :--- | :---: | :---: |
| **1** | **复数与矩阵代数** | 复数极坐标/欧拉公式、高斯消元求逆、完全主元数值秩 | 代数精确对齐 | PASS |
| **2** | **多项式代数与求根** | Horner 求值、卷积乘法、Durand-Kerner 复根全局迭代 | 代数精确对齐 | PASS |
| **3** | **一阶/二阶闭式解析解** | 欠阻尼 ($0<\zeta<1$)、临界 ($\zeta=1$)、过阻尼 ($\zeta>1$) 阶跃响应 | **$< 10^{-12} \ll 10^{-9}$** | PASS |
| **4** | **时域与积分性能指标** | 理论 $M_p, t_p, t_r, t_s$ 公式与轨迹实测 $\text{IAE}/\text{ISE}/\text{ITAE}$ | **$< 10^{-7} < 10^{-6}$** | PASS |
| **5** | **稳态误差与劳斯判据** | 终值定理 $K_p, K_v, K_a$ 与 Routh 阵列首列零元/全零行退化处理 | 判定一致率 100% | PASS |
| **6** | **定量频域稳定裕度** | $L(s)=\frac{2}{(s+1)^3}$ 二分法求解 $\omega_{pc}=\sqrt{3}, GM=12.04\text{dB}, PM=67.6^\circ$ | **误差 $< 10^{-3}$** | PASS |
| **7** | **状态空间与极点配置** | Kalman 能控/能观秩（识别零极点相消）与倒立摆 Ackermann 极点配置 | **误差 $< 10^{-9}$** | PASS |
| **8** | **无量纲收敛与后验误差证书** | 验证 $\mathcal{O}((|\lambda|dt)^4)$ 平均收敛阶 $4.0191$ 及 Richardson 后验误差证书 | **证书相对偏差 $< 0.25\%$** | PASS |
| **9** | **自适应步长控制 (RK45)** | 快极点高增益对象 $G(s)=\frac{100}{0.02s+1}$ 自动步长缩减与拒绝重算 | **局部误差 $< 5 \times 10^{-7}$** | PASS |
| **10** | **二维复射线与采样闭环谱半径** | 复射线稳定半径 $r_{\text{crit}}(\theta)$ 及采样闭环转移矩阵 $\rho(\mathbf{\Phi}_{\text{CL}}(dt)) < 1$ | **彻底消除稳定性反例** | PASS |
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
