# MoonControl 对标分析、应用场景与数值精度/稳定性边界说明

本文档参考 Python [`python-control`](https://python-control.readthedocs.io/)（基于 `SciPy`）、Julia [`ControlSystems.jl`](https://juliacontrol.github.io/ControlSystems.jl/) 等其他语言成熟控制系统库的设计范式，明确界定 `MoonControl` 的应用场景、与外部成熟库及 MoonBit 生态已有库（[`CAIMEOX/symbit/symphysics/control`](https://mooncakes.io/docs/CAIMEOX/symbit/symphysics/control)）的功能对比，以及**可信的数值精度与稳定性适用边界**。

---

## 一、典型应用场景

1. **《自动控制原理》与《现代控制工程》Web 交互式教学实验**
   - 依托 MoonBit 原生编译至 WebAssembly/JavaScript 的零依赖特性，配合仓库内置的 `web/index.html` 单文件实验台，学生与教师无需安装体积庞大的 MATLAB/Simulink 或配置 Python/SciPy 环境，即可在浏览器中实时开展二阶系统阻尼特性分析、劳斯表稳定性校验、Bode 裕度分析及 PID 参数整定实验。
2. **嵌入式与边缘控制回路离线原型仿真与参数预整定**
   - 面向直流伺服电机定位、温控加热器、DC-DC 电源环路等实际工程场景（普遍存在 PWM 占空比/电压饱和限幅、执行机构变化率限制、传感器/通信采样纯滞后 $e^{-\tau_d s}$ 以及数字控制器离散采样周期 $T_s$ 约束），在将控制算法烧录至 MCU 前，评估执行器饱和引起的积分 Windup 超调、整定反计算抗饱和增益 $K_b$，并校核采样周期 $T_s$ 是否引发数字环路振荡。
3. **MoonBit 原生仿真与机器人动力学基础组件**
   - 为纯 MoonBit 编写的轻量级物理仿真、无人机姿态控制环或 AGV 路径跟踪提供无 C FFI 依赖的 SISO 传递函数、可控标准型状态空间模型、状态反馈极点配置与 RK4/RK45 时域积分器。

---

## 二、与主流语言成熟库及 MoonBit 生态库的功能与边界对标

`MoonControl` 在 API 定位上参考了 Python `python-control` 与 `scipy.signal`/`scipy.integrate` 的核心子集，专注于**低阶单输入单输出（SISO）经典与现代控制分析及含常见执行器非线性的时域闭环仿真**。

| 功能模块 | Python `python-control` + `SciPy` | Julia `ControlSystems.jl` | MoonBit `CAIMEOX/symbit/symphysics/control` | MoonControl (`lycvvt/mooncontrol`) |
| :--- | :--- | :--- | :--- | :--- |
| **核心定位** | 工业/科研全功能数值控制库（底层依赖 Fortran SLICOT / LAPACK） | 高性能科学计算控制库（依赖 Julia LinearAlgebra / DifferentialEquations） | 符号计算控制模块（对齐 `sympy.physics.control`，基于 `@symcore.Expr`） | **纯 MoonBit 零 FFI 轻量级 SISO 数值控制与闭环仿真库** |
| **模型表示与代数** | SISO/MIMO `TransferFunction`, `StateSpace`, `FRD` | SISO/MIMO `TransferFunction`, `StateSpace`, `DelayLtiSystem` | SISO/MIMO 符号 `TransferFunction`, `StateSpace`, `Series`/`Parallel`/`Feedback` | SISO 数值 `TransferFunction`、可控标准型 `StateSpaceModel`、串/并/反馈连接 |
| **线性时域响应机制** | `scipy.linalg.expm` 矩阵指数离散递推 (`lsim` / `step_response`) | 矩阵指数精确 ZOH 离散化或 `OrdinaryDiffEq` | 符号部分分式展开 $\mathcal{L}^{-1}\{G(s)/s^k\}$ (`control_plots.mbt`) | 1/2 阶闭式解析解 + `matrix_exponential` 缩放平方矩阵指数 + RK4/RK45 数值积分 |
| **受限非线性闭环仿真** | 支持（通过 `control.input_output_response` + `solve_ivp`） | 支持（通过 `Simulink` 风格组件或 `DifferentialEquations.jl`） | **不支持**（纯线性拉氏变换，无法表示饱和限幅与积分抗饱和） | **支持 (`simulate_closed_loop_nonlinear`)**：执行器幅值限幅、速率限制、反算抗饱和 $K_b$ |
| **纯滞后环节 ($e^{-\tau_d s}$)** | Padé 有理逼近 (`control.pade`) 或时域缓冲 | 原生支持 `delay(tau)` 状态增广与缓冲仿真 | **不支持**（`control_plots.mbt:182` 抛出 `NotImplemented` 异常） | **支持**：时域 FIFO 环形延迟缓冲队列 + 阶跃负载扰动注入 |
| **频域稳定裕度** | `control.margin` / `stability_margins`（插值与数值求根） | `margin(sys)`（精确求解 $GM, PM, \omega_{gc}, \omega_{pc}$） | 仅输出离散频点幅相数组，无裕度求解函数 | **支持 (`compute_stability_margins`)**：复共轭相位解缠 + 50 步二分求根提取 $GM, PM, \omega_{gc}, \omega_{pc}, M_s$ |
| **能控能观与极点配置** | `ctrb`, `obsv`, `acker`, `place`（基于 SVD 与 Schur 分解） | `controllability`, `observability`, `place` | **不支持** | **支持 (`state_space.mbt`)**：全主元高斯消元秩判定 + `ackermann_pole_placement`（限低阶 $n \le 6$） |
| **大规模/MIMO/鲁棒控制** | 支持（$H_2/H_\infty$, LQR, Kalman, 阶数可达数百阶） | 支持（LQR, LQG, $H_\infty$, MPC） | 仅支持符号矩阵代数，无数值鲁棒综合 | **不支持（明确非目标）**：仅面向低阶 SISO 系统（建议阶数 $n \le 6$） |

---

## 三、可信的数值精度与稳定性边界（对标成熟数值库规范）

参考 `scipy.integrate.solve_ivp`（Hairer-Nørsett-Wanner 规范）与 `python-control` 的数值分析边界，本库将各算法的精度与稳定性保证严格区分为**确定性解析/矩阵级数边界**与**渐近启发式工程估计边界**，不作超越数值方法数学本质的过度承诺：

### 1. 确定性精度边界（仅限低阶理想 LTI 解析解与矩阵指数）
- **一阶与二阶标准系统闭式解析解 (`step_response_first_order`, `step_response_second_order`)**：
  - 直接由初等超越函数（$\exp, \sin, \cos, \sinh, \cosh$）闭式求值，不存在时间步进截断误差，在非奇异参数下（阻尼比 $|\zeta - 1| > 10^{-6}$ 或精确等于 $1$）精度受限于 IEEE-754 双精度浮点运算（实测误差 $< 10^{-12}$）。
- **小规模状态空间矩阵指数 (`matrix_exponential`)**：
  - 对标 `scipy.linalg.expm` 的缩放平方算法（Scaling and Squaring）：先选取 $2^s$ 将 $\|\mathbf{A} dt\|_\infty / 2^s$ 缩放至 $\le 0.5$，再采用 12 阶泰勒级数展开（局部级数截断余项 $\le \frac{0.5^{13}}{13!} < 2 \times 10^{-14}$），最后执行 $s$ 次平方还原。适用于低阶良态系统（$n \le 6$）的精确零阶保持（ZOH）线性离散步进。

### 2. ODE 数值积分的渐近误差估计边界（对标 `scipy.integrate.solve_ivp`）
对于含非线性饱和、速率限制或任意外部输入的常微分方程初值问题，本库提供定步长 RK4（`simulate_step`）与基于步长折半 Richardson 外推的自适应积分（`simulate_step_certified` / `simulate_step_adaptive`）：
- **渐近误差指示器（Asymptotic Error Indicator）而非严格区间上界**：
  - 步长折半后验估计式 $\hat{E} = \frac{16}{15}|y^{(h/2)} - y^{(h)}|$ 与 `solve_ivp(method='RK45')` 的嵌入式误差估计同属**渐近主项近似**，忽略了 $\mathcal{O}((\|\mathbf{A}\|h)^6)$ 高阶余项。
  - **有效前提条件**：
    1. **无量纲步长处于渐近收敛区**：必须满足 $\eta = |\lambda_{\max}| \cdot h \le 0.2$（在基准二阶系统 $G_A(s)=\frac{10}{s^2+2s+10}$ 上，当 $\eta \in [0.0158, 0.2530]$ 时实测平均收敛阶为 $4.0191$，渐近估计与真值误差相对偏差 $< 1\%$）；若固定秒级步长 $h$ 应用于快极点/高增益系统（如 $G(s)=\frac{100}{0.02s+1}$ 在 $h=0.01\text{s}$ 时 $\eta = 0.5$），单步误差可达 $2.40 \times 10^{-2}$，此时必须启用自适应步长缩减 $h$ 使其重新进入 $\eta \le 0.2$ 的渐近区。
    2. **右端函数局部光滑性 ($C^5$)**：在执行器刚刚触发饱和限幅（$\text{clamp}(u)$ 折点）或死区纯滞后阶跃跳变的单个步长内部，右端导数不连续，局部收敛阶会退化至 $\mathcal{O}(h^2) \sim \mathcal{O}(h)$，此时渐近误差估计仅作步长缩放启发式参考，不构成严格上界。
    3. **非刚性系统限制（Non-Stiff ODEs Only）**：显式 RK4/RK45 仅适用于刚度比 $S = \frac{\max |\text{Re}(\lambda_i)|}{\min |\text{Re}(\lambda_i)|} \le 100$ 的非刚性系统。当 $S > 100$ 时，步长受显式稳定域限制而非精度限制，应使用隐式积分器（本库 `check_system_boundaries` 将标记 `is_stiff = true` 告警）。

### 3. 稳定性判据与求根算法的数值条件数边界
- **连续系统劳斯判据 (`is_stable_routh`) 与极点求根 (`Polynomial::roots`)**：
  - 基于多项式系数的 Durand-Kerner 复根迭代与 Routh 表递推，受控于经典的 **Wilkinson 多项式条件数**。
  - **可信适用范围**：阶数 $n \le 6$、无高重数簇集根（重数 $\le 2$），且极点实部偏离虚轴 $|\text{Re}(\lambda_i)| > 10^{-6}$ 的良态系统。
  - **退化边界**：当极点极其接近虚轴（$|\text{Re}(\lambda_i)| \le 10^{-6}$，临界稳定）或阶数较高（$n > 8$）时，浮点舍入误差可能导致特征根跨过虚轴或劳斯表首列符号翻转，此时结果仅作临界预警，不能替代基于 QR/Schur 正交变换的工业级特征值求解器。
- **开环 RK4 方向稳定半径 (`rk4_directional_stability_radius`) 与采样闭环转移矩阵谱半径 (`closed_loop_pid_spectral_radius`)**：
  - 开环线性模态在 RK4 下收敛的条件为 $\text{Re}(\lambda_i) < 0$ 且 $|R(\lambda_i dt)| < 1$（沿复极点方向二分求解，最小稳定半径在 $\zeta \approx 0.75$ 处约为 $2.5657$，负实轴上约为 $2.7853$；若含右半平面极点则 `max_stable_dt = 0.0`）。
  - 对于数字采样 PID 闭环，零阶保持（ZOH）延迟与后向差分微分项会改变闭环特征根（例如 $G(s)=\frac{1}{s+1}$ 在 $K_p=8, K_i=4, K_d=0.2, dt=0.5\text{s}$ 时开环稳定但闭环发散 $\rho \approx 3.35 > 1$；当未滤波微分满足 $\mathbf{C}\mathbf{B}K_d > 1$ 时甚至对任意小步长 $dt \to 0^+$ 均存在模大于 1 的寄生极点）。
  - `closed_loop_pid_spectral_radius` 通过 Faddeev-LeVerrier 迹递推构造 $(n+2)$ 阶采样闭环转移矩阵 $\mathbf{\Phi}_{\text{CL}}(dt)$ 的特征多项式 $\det(z\mathbf{I} - \mathbf{\Phi}_{\text{CL}})$ 并求根估算 $\rho(\mathbf{\Phi}_{\text{CL}}(dt)) = \max_i |z_i|$。**其可信边界同样为低阶线性无饱和闭环（被控对象阶数 $n \le 6$）且特征根远离单位圆临界带（$||z_i| - 1| > 10^{-6}$）**；当存在执行器饱和非线性时，线性转移矩阵仅反映原点小信号局部线性化稳定性，不代表全局大信号渐近稳定。

---

## 四、基准算例实测汇总

通过执行 `moon test`（36 项单元测试）与 `moon run cmd/main` 可直接复现以下基准测试数据：

1. **二阶欠阻尼伺服系统 $G_A(s) = \frac{10}{s^2 + 2s + 10}$**：
   - 开环极点 $s = -1 \pm 3j$，超调量 $M_p = 35.09\%$；
   - 在执行器限幅 $u \in [-1.5, 1.5]$、纯滞后 $\tau_d = 0.02\text{s}$ 与 $t=2.5\text{s}$ 阶跃扰动 $d = -0.4$ 下，启用反算抗饱和（$K_b = 5.0$）将闭环超调由无抗饱和时的 $20.51\%$ 降低至 $5.81\%$，积分绝对误差 $\text{IAE}$ 由 $0.5842$ 降至 $0.5126$。
2. **倒立摆二阶不稳定对象 $G_B(s) = \frac{1}{s^2 - 9.8}$**：
   - 开环极点 $s = \pm 3.1305\text{ rad/s}$，Kalman 能控矩阵满秩（$\text{rank}(\mathcal{C}) = 2$）；通过 Ackermann 公式配置期望闭环极点 $s = -4 \pm 4j$，求得状态反馈增益 $\mathbf{K} = [41.8, 8.0]$。
3. **三阶滞后系统 $L(s) = \frac{2}{(s+1)^3}$ 频域裕度**：
   - 二分求根解得相位穿越频率 $\omega_{pc} = 1.73205\text{ rad/s}$（解析值 $\sqrt{3}$）、幅值裕度 $GM = 12.0412\text{ dB}$（解析值 $20\log_{10} 4$）、幅值穿越频率 $\omega_{gc} = 0.76642\text{ rad/s}$、相角裕度 $PM = 67.60^\circ$，灵敏度峰值 $M_s = 1.6667$。
