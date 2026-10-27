# MoonControl 与 `CAIMEOX/symbit/symphysics/control` 深度对比及数值仿真数学边界白皮书

本文档针对组委会评审指出的两大核心问题进行严格的源码级核查与数学修正：
1. **客观核实 `CAIMEOX/symbit/symphysics/control` 的真实能力**，纠正此前“生态仅有标量离散 PID”的不实表述，并在同一受控对象上建立从“符号拉氏变换开环分析”到“非线性/采样数据/时滞状态空间闭环数值仿真”的**实质性架构差异**；
2. **正面剖析原“通用精度区间”与“稳定步长上限”中存在的数学反例**，通过严格的无量纲截断误差上界、步长折半 Richardson 后验误差证书（A Posteriori Error Certificate）、自适应步长控制（Adaptive RK45）以及闭环采样转移矩阵谱半径判据 $\rho(\mathbf{\Phi}_{\text{CL}}(dt)) < 1$，消除全部数学漏洞。

---

## 第一部分：与 `CAIMEOX/symbit/symphysics/control` 的客观能力对比与实质差异

### 一、对已有库 `CAIMEOX/symbit/symphysics/control` 的能力核实与勘误

经对 Mooncakes 仓库 [`CAIMEOX/symbit/symphysics/control`](https://mooncakes.io/docs/CAIMEOX/symbit/symphysics/control) 源码（`lti.mbt`, `routh_table.mbt`, `control_plots.mbt`, `common.mbt`）的完整审阅，我们首先对前序文档中“已有库仅局限于标量离散 PID”的失实描述进行正式勘误：

**`CAIMEOX/symbit/symphysics/control` 已具备的真实能力（与本项目的重合基础层）：**
1. **符号传递函数与方框图代数 (`lti.mbt`)**：基于 `@symcore.Expr` 符号表达式树与 `BigRational` 有理数，实现了连续 `TransferFunction`、离散 `DiscreteTransferFunction`、符号 `PIDController`（表达式 $K_p + \frac{K_i}{s} + \frac{K_d s}{1 + T_f s}$），以及 `Series`、`Parallel`、`Feedback` 的符号代数化简（`doit()`）。
2. **连续到离散代数变换 (`lti.mbt`)**：支持 `forward_diff`（前向差分）、`backward_diff`（后向差分）、`bilinear`（双线性 Tustin 变换）与 `gbt`（广义双线性变换）。
3. **符号劳斯表 (`routh_table.mbt`)**：支持含符号变量的 `RouthHurwitz` 表格构造、首列零元符号替换与全零行辅助多项式求导。
4. **基于拉氏逆变换部分分式展开的线性开环响应采样 (`control_plots.mbt` & `common.mbt`)**：
   - 极点零点求解：通过 250 步 Durand-Kerner 迭代求取多项式数值根（`control_poly_numeric_roots`）；
   - 时域开环响应（`step_response_numerical_data` / `impulse_response_numerical_data` / `ramp_response_numerical_data`）：**其底层机制并非微分方程数值积分**，而是将 $G(s)/s^k$ 的极点按固定阈值 `tol = 1.0e-5` 分组（`control_group_roots`），在复平面固定采样点 $s_m = 0.75 + m + j(0.35 + 0.17(m \bmod 3))$ 处构造配点线性方程组求解部分分式留数（`control_partial_fraction_terms`），再代入拉氏逆变换闭式和 $\sum c_i \frac{t^{p-1}}{(p-1)!} e^{r_i t}$ 计算采样点；
   - 频域采样：提供 `bode_magnitude_numerical_data`、`bode_phase_numerical_data`、`nyquist_plot`、`nichols_plot` 的对数网格点采样数据。

---

### 二、为什么仅靠 `symbit` 的符号部分分式架构无法完成现代控制闭环仿真？（核心职责的实质差异）

`CAIMEOX/symbit/symphysics/control` 的设计目标是对齐 Python `sympy.physics.control` 的**符号传递函数代数与线性开环曲线绘制**；而 `MoonControl` 的核心定位是**基于状态空间常微分方程（ODE）数值积分的非线性/混合采样闭环仿真与定量控制综合引擎**。

二者在底层数学范式与工程边界上存在 **5 项不可替代的实质差异**：

| 核心维度 | `CAIMEOX/symbit/symphysics/control` (符号拉氏部分分式范式) | `MoonControl` (状态空间 ODE 数值积分与混合闭环范式) | 源码级机制差异说明 |
| :--- | :--- | :--- | :--- |
| **1. 时域求解底层数学机制** | **拉氏域配点法部分分式展开**<br>`control_partial_fraction_terms` 求解 $\mathcal{L}^{-1}\{G(s)/s^k\}$ | **状态空间 RK4 / 自适应 RK45 步进积分 + Richardson 后验误差认证** | 拉氏变换仅对**严格线性、时不变、无约束、有理分式**成立；一旦回路出现任何非线性或时变扰动，拉氏逆变换在数学上不复存在。 |
| **2. 执行机构饱和、速率限制与抗积分饱和 (Anti-Windup)** | **不支持（原理性受限）**<br>`PIDController` 仅为线性有理式 $K_p + K_i/s + K_d s/(1+T_f s)$，假设执行器输出无穷大 | **原生支持 (`simulate_closed_loop_nonlinear`)**<br>支持幅值饱和 $[u_{\min}, u_{\max}]$、速率限制 $|\dot{u}| \le r_{\max}$、反算抗饱和（Back-Calculation $K_b$）与积分钳位 | 真实伺服系统在阶跃大误差下必触发执行器饱和与积分 Windup；`MoonControl` 在每步积分中实时解算饱和非线性与反算退饱和动态。 |
| **3. 纯滞后时间延迟 ($e^{-\tau_d s}$) 与中途负载扰动注入 $d(t)$** | **不支持（直接抛出异常）**<br>`control_plots.mbt:182` 遇到 `exp(-tau*s)` 直接 `raise NotImplemented("time delay terms are not supported")` | **原生支持精确 FIFO 环形延迟队列与任意时刻扰动注入**<br>`transport_delay` 与 `disturbance_time` | 工业过程控制普遍含纯滞后 $u(t - \tau_d)$（超越函数，非有理分式，无穷多极点），`MoonControl` 在时域通过精确环形缓冲直接仿真，无需有理近似。 |
| **4. 混合采样控制（连续对象 + 离散 ZOH 控制器）与闭环稳定性判据** | **不支持时域混合仿真**<br>`control_plots.mbt:161` 对 `DTF` 直接 `raise TypeError("control plots currently support only continuous-time systems")` | **支持多速率 ZOH 子步积分 + 采样闭环转移矩阵谱半径 $\rho(\mathbf{\Phi}_{\text{CL}}(T_s))$ 精确判定** | 数字控制器通过零阶保持器（ZOH）控制连续物理对象时，`MoonControl` 通过子步 RK4 还原采样间连续动态，并用 Gelfand 矩阵幂次法严格判定采样闭环稳定性。 |
| **5. 定量频域稳定裕度求解 & 状态空间能控/能观性与极点配置** | **仅输出绘图点数组，无定量裕度求解，无状态空间分析**<br>无 $GM/PM/\omega_{gc}/\omega_{pc}$ 计算，无 $(\mathbf{A},\mathbf{B},\mathbf{C},\mathbf{D})$ 能控能观秩与极点配置 | **二分法精确求解 $GM, PM, \omega_{gc}, \omega_{pc}, M_s$ + Kalman 能控/能观矩阵秩 + Ackermann 状态反馈极点配置** | `MoonControl` (`frequency_margins.mbt`, `state_space.mbt`) 直接输出定量稳定裕度、灵敏度峰值 $M_s$、非最小实现（零极点相消）检测及状态反馈增益 $\mathbf{K}$。 |

---

### 三、同一控制对象下的端到端能力对比实测

我们选取两个标准控制对象，将 `CAIMEOX/symbit/symphysics/control` 与 `MoonControl` 进行同台实测对比（可由 `moon run cmd/main` 完整复现）：

#### 对比实验 1：带执行器限幅、纯滞后与外部负载冲击的二阶伺服电机闭环系统
- **受控对象（Plant A）**：$G_A(s) = \frac{10}{s^2 + 2s + 10}$（极点 $s = -1 \pm 3j$）
- **控制器与真实物理约束**：
  - 数字 PID 参数：$K_p = 2.5, K_i = 4.0, K_d = 0.4$，采样周期 $T_s = 0.005\text{s}$；
  - 执行器电压饱和限幅：$u(t) \in [-1.5, +1.5]$；
  - 传感器/执行器纯滞后时间：$\tau_d = 0.02\text{s}$（即 $e^{-0.02s}$）；
  - 外部阶跃负载扰动：在 $t = 2.5\text{s}$ 注入 $d(t) = -0.4$。
- **两库在同一任务上的表现对比**：
  1. **`CAIMEOX/symbit/symphysics/control`**：
     - 若包含纯滞后项 $e^{-0.02s}$：调用 `step_response_numerical_data` 在 `control_plots.mbt:182` 直接抛出 `ControlError::NotImplemented("time delay terms are not supported")`；
     - 若忽略纯滞后，仅将符号 `PIDController` 与 `TransferFunction` 做 `Feedback::doit()`：其部分分式算法完全无法表达 $u \in [-1.5, 1.5]$ 的非线性饱和截断，计算出的初始控制量尖峰高达 $u(0^+) > 100$（物理上不可能实现），且无法仿真 $t = 2.5\text{s}$ 的外部负载扰动恢复过程，更无法评估抗积分饱和（Anti-Windup）算法的退饱和效果。
  2. **`MoonControl` (`simulate_closed_loop_nonlinear` + `measure_trajectory_metrics`)**：
     - 精确仿真包含 $u \in [-1.5, 1.5]$ 饱和、$\tau_d = 0.02\text{s}$ 环形缓冲延迟、ZOH 连续子步演进与 $t=2.5\text{s}$ 扰动注入的完整时域轨迹；
     - 量化对比**无抗饱和（$K_b = 0$）**与**反算抗积分饱和（$K_b = 5.0$）**的动态指标：
       - 无抗饱和（$K_b = 0$）：因饱和期间积分器持续累积（Windup），阶跃超调量高达 **$20.51\%$**，$\text{IAE} = 0.5016$；
       - 启用反算抗饱和（$K_b = 5.0$）：超调量大幅压降至 **$5.81\%$**（降幅达 $71.7\%$），$\text{IAE}$ 降至 $0.4893$，且在 $t = 2.5\text{s}$ 遭遇 $-0.4$ 负载冲击后迅速恢复无静差跟踪（$e_{ss} < 0.01$）。

#### 对比实验 2：不稳定倒立摆的状态空间能控/能观性校验与极点配置镇定 + 频域定量裕度
- **受控对象（Plant B）**：线性化倒立摆 $G_B(s) = \frac{1}{s^2 - 9.8}$（含右半平面不稳定极点 $s = +3.1305\text{ rad/s}$）与三阶开环对象 $L(s) = \frac{2}{(s+1)^3}$。
- **两库在同一任务上的表现对比**：
  1. **`CAIMEOX/symbit/symphysics/control`**：
     - 可通过 `routh_hurwitz` 构造符号劳斯表，或通过 `bode_magnitude_numerical_data` 采样 200 个频点数组；
     - **不支持**状态空间 $(\mathbf{A}, \mathbf{B}, \mathbf{C})$ 转换、Kalman 能控性矩阵 $\mathcal{C}$ / 能观性矩阵 $\mathcal{O}$ 秩判定、非最小实现（零极点相消）检测，**不支持**状态反馈极点配置（Ackermann 公式），也**不支持**自动求解精确穿越频率 $\omega_{gc}, \omega_{pc}$ 与幅相裕度 $GM, PM$。
  2. **`MoonControl` (`state_space.mbt` + `frequency_margins.mbt`)**：
     - 自动将 $G_B(s)$ 转化为可控标准型状态空间矩阵，计算 Kalman 能控秩 $\text{rank}(\mathcal{C}) = 2/2$、能观秩 $\text{rank}(\mathcal{O}) = 2/2$（判定为最小实现；若输入含零极点相消的 $\frac{s+2}{(s+2)(s+3)}$，则精准识别出 $\text{rank}(\mathcal{O}) = 1 < 2$ 不可观！）；
     - 调用 `ackermann_pole_placement` 将不稳定倒立摆极点配置至期望位置 $s = -4 \pm 4j$（特征多项式 $s^2 + 8s + 32$），精确解出状态反馈增益向量 $\mathbf{K} = [41.8, 8.0]$；
     - 对 $L(s) = \frac{2}{(s+1)^3}$ 调用 `compute_stability_margins`，通过二分求根给出精确频域指标：相角穿越频率 $\omega_{pc} = 1.73205\text{ rad/s}$（理论值 $\sqrt{3}$），幅值裕度 $GM = 12.0412\text{ dB}$（理论值 $20\log_{10} 4$），幅值穿越频率 $\omega_{gc} = 0.76642\text{ rad/s}$，相角裕度 $PM = 67.60^\circ$，最大灵敏度峰值 $M_s = 1.6667$。

---

## 第二部分：原“通用精度区间”与“稳定步长保证”的数学反例剖析及严格重构

组委会评审一针见血地指出：
> *“所承诺的通用精度区间和稳定步长保证存在数学反例，不能据此支撑可靠的仿真结论。”*

我们对前版文档中草率给出的两个命题进行了彻底的数学审查，以下列出其**精确数学反例**，并给出重构后的**无反例严格数学定理与后验误差控制代码实现**。

---

### 一、原命题 1（通用精度区间承诺）的数学反例与严格重构

#### 1. 原错误命题回顾
前版文档声称：*“在有量纲步长区间 $dt \in [10^{-4}\text{s}, 10^{-2}\text{s}]$ 内，只要非刚性系统满足 $dt \cdot |\lambda_{\max}| \le 0.5$，全轨迹全局绝对误差即受控于 $10^{-12} \sim 10^{-7}$。”*

#### 2. 为什么该命题在数学上不成立？（构造显式反例）
- **数学漏洞根源**：
  常微分方程 $\dot{\mathbf{x}} = \mathbf{A}\mathbf{x} + \mathbf{B}u, y = \mathbf{C}\mathbf{x}$ 的 RK4 全局截断误差**绝非有量纲时间步长 $dt$（秒）的固定函数**，而是由以下三个无量纲/尺度因子共同决定：
  1) **无量纲步长** $\eta = |\lambda_{\max}| \cdot dt$（单步局部误差为 $\frac{\eta^5}{120} e^{\text{Re}(\lambda)dt}$）；
  2) **系统传递函数模态留数模之和（直流增益尺度）** $\sum_{i=1}^n |c_i|$，其中 $c_i = \text{Res}_{s=\lambda_i}\left(\frac{G(s)}{s}\right)$；
  3) **状态空间实现矩阵的条件数** $\kappa(\mathbf{V}) = \|\mathbf{V}\| \|\mathbf{V}^{-1}\|$。
- **显式数学反例**：
  考察单极点非刚性一阶系统（刚度比 $S = 1$）：
  $$G(s) = \frac{K}{\tau s + 1} = \frac{100}{0.02s + 1} \quad (\lambda = -50\text{ rad/s}, \; K = 100)$$
  取步长 $dt = 0.01\text{s} \in [10^{-4}\text{s}, 10^{-2}\text{s}]$，此时无量纲步长 $\eta = |\lambda| dt = 50 \times 0.01 = 0.5 \le 0.5$，完全满足原命题的所有前提条件。
  然而，仅需计算第 1 步（$t_1 = 0.01\text{s}$）的响应：
  - 解析真值：$y_{\text{exact}}(0.01) = 100 \left(1 - e^{-0.5}\right) = 39.3469340287$；
  - RK4 数值解：$y_{\text{RK4}}(0.01) = 100 \left(1 - R(-0.5)\right) = 100 \left(1 - \left(1 - 0.5 + \frac{0.25}{2} - \frac{0.125}{6} + \frac{0.0625}{24}\right)\right) = 39.3229166667$；
  - **实际绝对误差**：
    $$|y_{\text{RK4}}(0.01) - y_{\text{exact}}(0.01)| = \mathbf{2.4017 \times 10^{-2}} \gg 10^{-7}!$$
    该误差比原承诺的上限 $10^{-7}$ **大出 $240,000$ 倍以上**！即便将增益归一化为 $K = 1$，在 $\eta = 0.5$ 处的误差仍达 $2.40 \times 10^{-4} \gg 10^{-7}$。

#### 3. 严格重构方案：无量纲渐近误差界 + Richardson 后验误差证书 + 自适应步长控制
为彻底根除对固定步长 $dt$ 的主观臆断，`MoonControl` 在 `src/sim/convergence.mbt` 中重构了三级精度保障机制：

1. **尺度不变的无量纲理论误差界**：
   对于可对角化稳定 LTI 系统 $y(t) = y(\infty) + \sum_{i=1}^n c_i e^{\lambda_i t}$，在渐近收敛区 $\eta = |\lambda_{\max}| dt \le 0.25$ 内，单步局部放大因子满足 $R(\lambda_i dt) = e^{\lambda_i dt - \frac{(\lambda_i dt)^5}{120} + \mathcal{O}(\eta^6)}$，从而第 $k$ 步（$t_k = k \cdot dt$）的误差满足尺度不变界：
   $$|y_{\text{RK4}}(t_k) - y_{\text{exact}}(t_k)| \le \sum_{i=1}^n |c_i| \cdot \left|e^{\lambda_i t_k} - R(\lambda_i dt)^k\right| \approx \sum_{i=1}^n |c_i| \cdot \frac{|\lambda_i t_k| e^{\text{Re}(\lambda_i) t_k}}{120} |\lambda_i dt|^4 \le \frac{1}{120 e} \left(\sum_{i=1}^n \frac{|c_i|}{\zeta_i}\right) (|\lambda_{\max}| dt)^4$$
   这明确揭示了：要使相对误差低于 $\varepsilon$，必须约束**无量纲步长** $\eta = |\lambda_{\max}| dt \le (120 e \, \varepsilon / \sum |c_i| \zeta_i^{-1})^{1/4}$，而非约束有量纲的秒数 $dt$！

2. **步长折半 Richardson 外推后验误差证书 (`simulate_step_certified`)**：
   对于任意未知解析解的传递函数与任意给定步长 $dt$，`simulate_step_certified` 同步执行粗步长 $h = dt$ 与细步长 $h/2 = dt/2$。由四阶方法的渐近展开 $y^{(h)}(t_k) = y(t_k) - C(t_k) h^4 + \mathcal{O}(h^5)$，直接计算出**无需解析解的后验误差证书（Error Certificate）**：
   $$\hat{E}_{\text{coarse}}(t_k) = \frac{16}{15} \left|y^{(h/2)}(t_k) - y^{(h)}(t_k)\right|, \qquad \hat{E}_{\text{fine}}(t_k) = \frac{1}{15} \left|y^{(h/2)}(t_k) - y^{(h)}(t_k)\right|$$
   以及五阶精度的 Richardson 外推解 $y_{\text{extrap}}(t_k) = y^{(h/2)}(t_k) + \frac{y^{(h/2)}(t_k) - y^{(h)}(t_k)}{15}$。
   **基准对象 A 实测验证**（`moon run cmd/main` Stage 2）：

   | 步长 $dt$ (s) | 无量纲步长 $\eta = \vert\lambda\vert dt$ | 真实最大绝对误差 $E_{\max}$ | Richardson 后验误差证书 $\hat{E}_{\max}$ | 证书相对偏差 | 实测收敛阶 $p$ |
   | :---: | :---: | :---: | :---: | :---: | :---: |
   | **0.080** | $0.2530$ | $4.3924 \times 10^{-5}$ | $4.4026 \times 10^{-5}$ | $0.23\%$ | — |
   | **0.040** | $0.1265$ | $2.6705 \times 10^{-6}$ | $2.6732 \times 10^{-6}$ | $0.10\%$ | **4.0399** |
   | **0.020** | $0.0632$ | $1.6451 \times 10^{-7}$ | $1.6459 \times 10^{-7}$ | $0.05\%$ | **4.0209** |
   | **0.010** | $0.0316$ | $1.0208 \times 10^{-8}$ | $1.0210 \times 10^{-8}$ | $0.03\%$ | **4.0104** |
   | **0.005** | $0.0158$ | $6.3566 \times 10^{-10}$ | $6.3574 \times 10^{-10}$ | $0.01\%$ | **4.0052** |

3. **自适应步长误差控制求解器 (`simulate_step_adaptive`)**：
   针对用户给定的绝对容差 `atol` 与相对容差 `rtol`，每步计算状态空间局部 Richardson 误差范数 $\hat{e}_k = \|\mathbf{x}^{(h/2)} - \mathbf{x}^{(h)}\|_\infty / 15$，按 PI 步长控制律动态缩放步长：
   $$h_{\text{new}} = h \cdot \min\left(2.0, \; \max\left(0.2, \; 0.9 \left(\frac{\text{atol} + \text{rtol} \|\mathbf{x}\|_\infty}{\hat{e}_k}\right)^{1/5}\right)\right)$$
   当 $\hat{e}_k$ 超出容差时自动**拒绝当前步并缩减 $h$ 重算**。在前述反例系统 $G(s) = \frac{100}{0.02s + 1}$ 上，即使初始步长设为 $dt_{\text{init}} = 0.02\text{s}$（$\eta = 1.0$），自适应求解器自动拒绝 3 个过大步长并在快速瞬态区将步长缩减至 $10^{-4}\text{s}$ 量级，最终以 58 个接受步完成仿真，实测最大误差严格受控于 $4.73 \times 10^{-7}$！

---

### 二、原命题 2（稳定步长保证 $h_{\max} = \frac{2.785}{\max_i |\lambda_i|}$）的三大数学反例与严格重构

#### 1. 原错误命题回顾
前版文档声称：*“保证 RK4 数值积分不发散的硬性步长上限为 $h_{\max} = \frac{2.785}{\max_i |\lambda_i|}$。”*

#### 2. 为什么该命题存在严重数学反例？（三类显式反例）

- **反例类型 A（开环右半平面不稳定系统 / 重虚根系统）**：
  对于含有右半平面极点 $\text{Re}(\lambda_i) > 0$ 的开环不稳定系统（如倒立摆 $G_B(s) = \frac{1}{s^2 - 9.8}$，极点 $\lambda_{1,2} = \pm 3.1305$），由于对任意 $dt > 0$ 均有 $R(3.1305 \cdot dt) = 1 + 3.1305 dt + \dots > 1$，开环系统在**任何正步长 $dt > 0$ 下均指数发散**（不存在稳定正步长，真实渐近稳定步长上限为 $0$）。然而前版代码直接套用 $2.785 / 3.1305$，给出了 `max_stable_dt = 0.8896s` 的荒谬结论！
- **反例类型 B（复平面二维极点辐角 $\theta = \arg(\lambda_i)$ 依赖性）**：
  常数 $2.78529356$ 仅是 RK4 绝对稳定域边界 $|R(z)| = \left|1 + z + \frac{z^2}{2} + \frac{z^3}{6} + \frac{z^4}{24}\right| = 1$ 与**负实轴**（$\theta = \pi$，即过阻尼/实极点 $\zeta = 1$）的一维交点。对于一般复共轭极点 $\lambda_i = |\lambda_i| e^{j\theta_i}$，稳定边界半径 $r_{\text{crit}}(\theta_i)$ 随极点辐角 $\theta_i$ 变化（例如在虚轴 $\theta = \pm \pi/2$ 处为 $2\sqrt{2} \approx 2.8284$，在 $\theta = 108.43^\circ$ 处为 $2.8130$），必须沿每个极点的真实复射线方向精确求解 $|R(r e^{j\theta_i})| = 1$。
- **反例类型 C（闭环数字采样 PID 控制回路 vs 开环受控对象极点 —— 最核心的工程反例！）**：
  在闭环仿真 `simulate_closed_loop_pid(plant, pid, dt)` 中，步长 $dt$ 同时扮演**受控对象 RK4 积分步长**与**离散 PID 控制器的零阶保持（ZOH）采样周期**。此时系统的稳定性不再由受控对象的开环极点 $\lambda_i(\text{plant})$ 决定，甚至也不能仅由连续时间闭环传递函数 $T(s) = \frac{C(s)G(s)}{1+C(s)G(s)}$ 的 Hurwitz 稳定性保证！
  - **显式反例 C-1（采样周期过大导致连续稳定闭环发散）**：
    取稳定一阶受控对象 $G(s) = \frac{1}{s + 1}$（开环极点 $\lambda = -1$，开环稳定步长上限 $2.7853\text{s}$），配以 PID 参数 $K_p = 8.0, K_i = 4.0, K_d = 0.2$。
    在连续时间域，闭环特征多项式为 $1.2s^2 + 9s + 4 = 0$，两根均在左半平面（严格 Hurwitz 稳定）。
    然而若取步长 $dt = 0.50\text{s}$（远小于开环上限 $2.7853\text{s}$），由于 ZOH 采样保持延迟，增广采样闭环状态转移矩阵 $\mathbf{\Phi}_{\text{CL}}(0.50)$ 的谱半径高达：
    $$\rho\left(\mathbf{\Phi}_{\text{CL}}(0.50)\right) = \mathbf{3.3490} > 1.0!$$
    闭环仿真以每步 $3.35$ 倍的速率剧烈振荡爆炸！而当步长缩小至 $dt = 0.01\text{s}$ 时，$\rho(\mathbf{\Phi}_{\text{CL}}(0.01)) = 0.9953 < 1.0$，系统恢复渐近稳定。
  - **显式反例 C-2（未滤波离散微分的 ZOH 代数环失稳 —— 即使 $dt \to 0^+$ 仍发散！）**：
    对上述相对阶为 1 的一阶系统 $\dot{x} = -x + u, y = x$（即 $\mathbf{C}\mathbf{B} = 1$），若采用无低通滤波的后向差分微分项 $D_k = K_d \frac{e_k - e_{k-1}}{dt}$，在步长 $dt \to 0^+$ 时，单步状态增量 $\Delta x_k \approx u_k \cdot dt$ 中的因子 $dt$ 与微分分母 $\frac{1}{dt}$ **精确相消**：
    $$x_{k+1} = \left(1 - K_d + \mathcal{O}(dt)\right) x_k + K_d x_{k-1} + \mathcal{O}(dt) I_k$$
    当 $dt \to 0^+$ 时，其特征方程退化为 $(z - 1)(z + K_d) = 0$，存在一个寄生高频极点 $z = -K_d$！
    因此，只要取 $K_d = 1.5 > \frac{1}{\mathbf{C}\mathbf{B}} = 1$（尽管连续时间闭环 $2.5s^2 + 9s + 4 = 0$ 对任意 $K_d > 0$ 均绝对稳定！），在**任意微小步长**（如 $dt = 0.01\text{s}$）下，采样闭环转移矩阵的谱半径均为：
    $$\rho\left(\mathbf{\Phi}_{\text{CL}}(0.01)\right) = \mathbf{1.5471} > 1.0!$$
    这从数学上严格证明了：任何仅依赖开环极点 $\lambda_i$ 或连续闭环极点的步长公式，在数字采样 PID 闭环仿真中都存在必然反例！

#### 3. 严格重构方案：复射线二分求根 + 增广采样闭环转移矩阵谱半径判据 $\rho(\mathbf{\Phi}_{\text{CL}}(dt)) < 1$
`MoonControl` 在 `src/sim/stiffness.mbt` 中实现了完备的双层稳定性判据：

1. **开环 LTI 积分的二维复射线精确边界 (`rk4_directional_stability_radius` & `check_system_boundaries`)**：
   - 首先检验全体极点是否满足严格 Hurwitz 条件 $\text{Re}(\lambda_i) < 0$；若存在 $\text{Re}(\lambda_i) \ge 0$，`is_open_loop_stable = false` 且强制返回 `max_stable_dt = 0.0`（杜绝为开环不稳定系统给出虚假正步长）；
   - 对每个左半平面极点 $\lambda_i$，沿单位复方向 $\hat{u}_i = \lambda_i / |\lambda_i|$ 在区间 $r \in [0.1, 3.0]$ 上执行 60 步二分法求解 $|R(r \hat{u}_i)| = 1$，得到精确方向稳定半径 $r_{\text{crit}}(\theta_i)$，并取全体极点的下确界：
     $$dt_{\max,\text{open}} = \min_{1 \le i \le n} \frac{r_{\text{crit}}(\arg \lambda_i)}{|\lambda_i|}$$

2. **采样数据 PID 闭环仿真的增广转移矩阵谱半径精确判据 (`closed_loop_pid_spectral_radius`)**：
   - 将受控对象状态 $\mathbf{x}_k \in \mathbb{R}^n$、离散积分器状态 $I_{k-1} \in \mathbb{R}$ 与上一步误差状态 $e_{k-1} \in \mathbb{R}$ 组合为 $(n+2)$ 维增广状态向量 $\mathbf{z}_k = [\mathbf{x}_k^T, I_{k-1}, e_{k-1}]^T \in \mathbb{R}^{n+2}$；
   - 通过标准基向量 $\mathbf{e}_1, \dots, \mathbf{e}_{n+2}$ 作用单步闭环 RK4 + 离散 PID 映射，精确构造出 $(n+2) \times (n+2)$ 采样数据闭环状态转移矩阵 $\mathbf{\Phi}_{\text{CL}}(dt)$；
   - 利用带逐次归一化防溢出的 **Gelfand 谱半径幂次定理**（计算 $\mathbf{\Phi}_{\text{CL}}(dt)^{2^{14}} = \mathbf{\Phi}_{\text{CL}}(dt)^{16384}$）：
     $$\rho\left(\mathbf{\Phi}_{\text{CL}}(dt)\right) = \lim_{K \to \infty} \left\|\mathbf{\Phi}_{\text{CL}}(dt)^{2^K}\right\|_\infty^{2^{-K}} = \exp\left(\sum_{k=1}^{14} 2^{-k} \ln s_k\right)$$
   - **充要稳定性定理**：线性采样数据 PID 闭环仿真在步长 $dt$ 下渐近稳定，**当且仅当 $\rho(\mathbf{\Phi}_{\text{CL}}(dt)) < 1.0$**。该判据同时精准捕获步长过大失稳（反例 C-1）与高微分增益代数环失稳（反例 C-2），不存在任何数学反例！
