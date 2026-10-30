# MoonControl 控制理论公式推导与算法笔记

本文档系统梳理了 **MoonControl** 核心算法的数学推导与理论基础，供开发者与赛事评审核对参考。

---

## 一、一阶系统单位阶跃响应推导

标准一阶惯性环节传递函数为：
$$G(s) = \frac{K}{Ts + 1}$$

其中 $K$ 为静态放大系数（直流增益），$T$ 为时间常数。

当输入为单位阶跃信号 $r(t) = 1(t) \iff R(s) = \frac{1}{s}$ 时，系统输出的拉普拉斯变换为：
$$Y(s) = G(s) R(s) = \frac{K}{s(Ts + 1)} = \frac{K/T}{s(s + 1/T)}$$

进行部分分式展开：
$$Y(s) = \frac{A}{s} + \frac{B}{s + 1/T}$$
其中：
$$A = \lim_{s \to 0} s Y(s) = K$$
$$B = \lim_{s \to -1/T} (s + 1/T) Y(s) = -K$$

对 $Y(s) = K \left( \frac{1}{s} - \frac{1}{s + 1/T} \right)$ 进行拉普拉斯反变换，得到时域解析解：
$$\mathbf{y(t) = K \left(1 - e^{-t/T}\right) \quad (t \ge 0)}$$

### 关键特征点验证：
- 当 $t = 0$ 时：$y(0) = K(1 - 1) = 0$
- 当 $t = T$ 时：$y(T) = K(1 - e^{-1}) \approx 0.63212 K$
- 当 $t = 2T$ 时：$y(2T) = K(1 - e^{-2}) \approx 0.86466 K$
- 当 $t = 3T$ 时：$y(3T) = K(1 - e^{-3}) \approx 0.95021 K$
- 当 $t \to \infty$ 时：$y(\infty) = K$

在 `src/sim/sim_wbtest.mbt` 中，数值仿真 RK4 的计算结果与上述解析公式的最大绝对误差在 $dt=0.0005$ 时为 $7.77 \times 10^{-15}$，并严格遵循四阶收敛律 $\mathcal{O}(dt^4)$。

---

## 二、标准二阶系统单位阶跃响应推导

典型二阶系统闭环传递函数为：
$$G(s) = \frac{\omega_n^2}{s^2 + 2\zeta\omega_n s + \omega_n^2}$$
其中 $\omega_n$ 为无阻尼固有角频率，$\zeta$ 为阻尼比。特征方程根为：
$$s_{1,2} = -\zeta\omega_n \pm \omega_n \sqrt{\zeta^2 - 1}$$

输出信号拉氏变换为：
$$Y(s) = \frac{1}{s} \cdot \frac{\omega_n^2}{s^2 + 2\zeta\omega_n s + \omega_n^2}$$

### 1. 欠阻尼情况 ($0 < \zeta < 1$)

特征根为一对共轭复数：
$$s_{1,2} = -\zeta\omega_n \pm j \omega_d, \quad \omega_d = \omega_n \sqrt{1 - \zeta^2}$$
输出表达式展开为：
$$Y(s) = \frac{1}{s} - \frac{s + 2\zeta\omega_n}{(s + \zeta\omega_n)^2 + \omega_d^2} = \frac{1}{s} - \frac{s + \zeta\omega_n}{(s + \zeta\omega_n)^2 + \omega_d^2} - \frac{\zeta\omega_n}{\omega_d} \frac{\omega_d}{(s + \zeta\omega_n)^2 + \omega_d^2}$$

拉普拉斯反变换：
$$y(t) = 1 - e^{-\zeta\omega_n t} \cos(\omega_d t) - \frac{\zeta}{\sqrt{1 - \zeta^2}} e^{-\zeta\omega_n t} \sin(\omega_d t)$$
引入辅助角 $\phi = \arccos(\zeta) = \arctan\left(\frac{\sqrt{1-\zeta^2}}{\zeta}\right)$：
$$\cos(\phi) = \zeta, \quad \sin(\phi) = \sqrt{1 - \zeta^2}$$
合并正余弦项得到标准解析解：
$$\mathbf{y(t) = 1 - \frac{e^{-\zeta\omega_n t}}{\sqrt{1 - \zeta^2}} \sin(\omega_d t + \phi)}$$

### 2. 临界阻尼情况 ($\zeta = 1$)

特征根为二重负实根 $s_{1,2} = -\omega_n$：
$$Y(s) = \frac{\omega_n^2}{s(s + \omega_n)^2} = \frac{1}{s} - \frac{1}{s + \omega_n} - \frac{\omega_n}{(s + \omega_n)^2}$$
拉氏反变换得到：
$$\mathbf{y(t) = 1 - e^{-\omega_n t}(1 + \omega_n t)}$$

### 3. 过阻尼情况 ($\zeta > 1$)

特征根为两个相异负实根 $p_1, p_2 < 0$：
$$p_{1,2} = -\omega_n(\zeta \mp \sqrt{\zeta^2 - 1})$$
展开为单极点：
$$Y(s) = \frac{1}{s} + \frac{A}{s - p_1} + \frac{B}{s - p_2}$$
留数求解：$A = \frac{p_2}{p_1 - p_2}, B = \frac{-p_1}{p_1 - p_2}$，故：
$$\mathbf{y(t) = 1 + \frac{p_2}{p_1 - p_2} e^{p_1 t} - \frac{p_1}{p_1 - p_2} e^{p_2 t}}$$

---

## 三、欠阻尼二阶系统动态性能指标推导

### 1. 峰值时间 $t_p$
令导数 $\frac{dy(t)}{dt} = 0$：
$$\frac{dy}{dt} = \frac{\omega_n}{\sqrt{1-\zeta^2}} e^{-\zeta\omega_n t} \sin(\omega_d t) = 0$$
首个非零正根满足 $\omega_d t_p = \pi$：
$$\mathbf{t_p = \frac{\pi}{\omega_d} = \frac{\pi}{\omega_n \sqrt{1 - \zeta^2}}}$$

### 2. 最大超调量 $M_p$
将 $t_p$ 代入 $y(t)$：
$$y(t_p) = 1 - \frac{e^{-\zeta\omega_n \pi / \omega_d}}{\sqrt{1-\zeta^2}} \sin(\pi + \phi) = 1 + e^{-\frac{\pi\zeta}{\sqrt{1-\zeta^2}}}$$
超调量定义为峰值超出稳态值的比例：
$$\mathbf{M_p = e^{-\frac{\pi\zeta}{\sqrt{1-\zeta^2}}}, \quad \sigma\% = M_p \times 100\%}$$

### 3. 上升时间 $t_r$
定义为响应曲线由 0 首次到达稳态值 1 的时刻：
$$y(t_r) = 1 \implies \sin(\omega_d t_r + \phi) = 0 \implies \omega_d t_r + \phi = \pi$$
$$\mathbf{t_r = \frac{\pi - \phi}{\omega_d}}$$

### 4. 调节时间 $t_s$
输出曲线进入并稳定在稳态值允许误差带（$\pm \Delta$）内所需的最短时间，由包络线衰减控制：
$$\frac{e^{-\zeta\omega_n t_s}}{\sqrt{1-\zeta^2}} \le \Delta$$
- 2% 准则 ($\Delta = 0.02$)：$e^{-\zeta\omega_n t_s} \approx 0.02 \implies \mathbf{t_s \approx \frac{4}{\zeta\omega_n}}$
- 5% 准则 ($\Delta = 0.05$)：$e^{-\zeta\omega_n t_s} \approx 0.05 \implies \mathbf{t_s \approx \frac{3}{\zeta\omega_n}}$

---

## 四、终值定理与稳态误差 $e_{ss}$

单位负反馈系统的误差传递函数为：
$$E(s) = \frac{R(s)}{1 + G(s)}$$
由拉普拉斯变换终值定理：
$$e_{ss} = \lim_{t \to \infty} e(t) = \lim_{s \to 0} s E(s) = \lim_{s \to 0} \frac{s R(s)}{1 + G(s)}$$

设开环传递函数为 $G(s) = \frac{K \prod(s - z_i)}{s^N \prod(s - p_j)}$，其中 $N$ 为系统型别：
- **0 型系统 ($N = 0$)**：
  - 阶跃输入 ($R(s) = 1/s$)：$e_{ss} = \frac{1}{1 + K_p}$，位置误差常数 $K_p = \lim_{s \to 0} G(s)$
  - 斜坡输入 ($R(s) = 1/s^2$)：$e_{ss} = \infty$
- **1 型系统 ($N = 1$)**：
  - 阶跃输入：$e_{ss} = 0$
  - 斜坡输入：$e_{ss} = \frac{1}{K_v}$，速度误差常数 $K_v = \lim_{s \to 0} s G(s)$
  - 抛物线输入：$e_{ss} = \infty$
- **2 型系统 ($N = 2$)**：
  - 阶跃、斜坡输入：$e_{ss} = 0$
  - 抛物线输入：$e_{ss} = \frac{1}{K_a}$，加速度误差常数 $K_a = \lim_{s \to 0} s^2 G(s)$

---

## 五、Routh-Hurwitz 代数判据与退化处理

对于特征多项式 $D(s) = a_n s^n + a_{n-1} s^{n-1} + \dots + a_0 = 0$：
Routh 阵列定义递推公式：
$$c_{i, j} = \frac{c_{i-1, 0} c_{i-2, j+1} - c_{i-2, 0} c_{i-1, j+1}}{c_{i-1, 0}}$$

### 边界退化处理算法：
1. **首列元素为 0 但所在行不全为 0**：
   引入微小正扰动 $\epsilon > 0$（代码中取 $\epsilon = 10^{-6}$），替换 $c_{i-1, 0}$ 继续推导后续行。
2. **整行为 0（全零行）**：
   出现于特征根具有关于原点对称（包括纯虚根 $\pm j\omega$ 或对称实根 $\pm \sigma$）时。
   构建上一行的辅助多项式 $A(s) = c_{i-2, 0} s^k + c_{i-2, 1} s^{k-2} + \dots$，对其求导 $\frac{dA(s)}{ds}$，用导数系数填补全零行。全零行表明系统临界稳定或不稳定。

系统严格稳定的充要条件：Routh 阵列第一列所有元素符号严格相同，且无全零行。

---

## 六、可控标准型状态空间实现与 RK4 积分

对于任意 $n$ 阶严格适定传递函数：
$$G(s) = \frac{b_{n-1}s^{n-1} + \dots + b_1 s + b_0}{s^n + a_{n-1}s^{n-1} + \dots + a_0}$$
选取状态变量 $x_1 = y, x_2 = \dot{x}_1, \dots, x_n = \dot{x}_{n-1}$，建立可控标准型微分方程组：
$$\begin{cases}
\dot{x}_1 = x_2 \\
\dot{x}_2 = x_3 \\
\vdots \\
\dot{x}_{n-1} = x_n \\
\dot{x}_n = -a_0 x_1 - a_1 x_2 - \dots - a_{n-1} x_n + u(t)
\end{cases}$$
输出方程：$y(t) = b_0 x_1 + b_1 x_2 + \dots + b_{n-1} x_n$。

经典四阶 Runge-Kutta (RK4) 每步演进公式：
$$k_1 = f(x_k, u)$$
$$k_2 = f(x_k + \frac{dt}{2} k_1, u)$$
$$k_3 = f(x_k + \frac{dt}{2} k_2, u)$$
$$k_4 = f(x_k + dt \cdot k_3, u)$$
$$x_{k+1} = x_k + \frac{dt}{6}(k_1 + 2 k_2 + 2 k_3 + k_4)$$

---

## 七、无量纲数值误差界、步长折半渐近误差估计与闭环采样稳定分析

### 1. 无量纲局部与全局截断误差推导
对于可对角化稳定 LTI 系统 $y(t) = y(\infty) + \sum_{i=1}^n c_i e^{\lambda_i t}$（其中 $c_i = \text{Res}_{s=\lambda_i}(G(s)/s)$），四阶龙格-库塔算法单步放大因子 $R(z) = 1 + z + \frac{z^2}{2} + \frac{z^3}{6} + \frac{z^4}{24}$ 满足 $R(\lambda_i h) = e^{\lambda_i h} - \frac{(\lambda_i h)^5}{120} + \mathcal{O}(|\lambda_i h|^6)$。
在渐近收敛区（$\eta = |\lambda_{\max}| h \le 0.2$）内，全局截断误差取决于**无量纲步长** $\eta = |\lambda_{\max}| h$ 与**模态留数模之和** $\sum |c_i|$，而非单纯有量纲秒数 $h$：
$$|y_{\text{RK4}}(t_k) - y_{\text{exact}}(t_k)| \le \sum_{i=1}^n |c_i| \cdot \left|e^{\lambda_i t_k} - R(\lambda_i h)^k\right| \lesssim \frac{1}{120 e} \left(\sum_{i=1}^n \frac{|c_i|}{\zeta_i}\right) (|\lambda_{\max}| h)^4$$

### 2. 步长折半 Richardson 渐近误差指示器与自适应步长控制
参考 `scipy.integrate.solve_ivp` 的渐近局部误差控制方法，`simulate_step_richardson` 同步运行步长 $h$ 与半步长 $h/2$，在渐近收敛区（$\eta \le 0.2$）且右端光滑的前提下，利用渐近展开 $y^{(h)}(t_k) = y(t_k) - C(t_k) h^4 + \mathcal{O}(h^5)$ 给出渐近后验误差估计（Asymptotic Error Indicator，用于工程步长缩放而非无条件区间上界）：
$$\hat{E}_{\text{coarse}}(t_k) = \frac{16}{15} \left|y^{(h/2)}(t_k) - y^{(h)}(t_k)\right|, \qquad \hat{E}_{\text{fine}}(t_k) = \frac{1}{15} \left|y^{(h/2)}(t_k) - y^{(h)}(t_k)\right|$$
在基准受控对象 $G(s) = \frac{10}{s^2 + 2s + 10}$ 上实测步长敏感度与 Richardson 渐近估计对比：
- $dt = 0.0800\text{s} \; (\eta = 0.2530)$: $E_{\max} = 4.3924 \times 10^{-5}$, $\hat{E}_{\max} = 4.4026 \times 10^{-5}$
- $dt = 0.0400\text{s} \; (\eta = 0.1265)$: $E_{\max} = 2.6705 \times 10^{-6}$, $\hat{E}_{\max} = 2.6732 \times 10^{-6} \implies p = 4.04$
- $dt = 0.0200\text{s} \; (\eta = 0.0632)$: $E_{\max} = 1.6451 \times 10^{-7}$, $\hat{E}_{\max} = 1.6459 \times 10^{-7} \implies p = 4.02$
- $dt = 0.0100\text{s} \; (\eta = 0.0316)$: $E_{\max} = 1.0208 \times 10^{-8}$, $\hat{E}_{\max} = 1.0210 \times 10^{-8} \implies p = 4.01$
- $dt = 0.0050\text{s} \; (\eta = 0.0158)$: $E_{\max} = 6.3566 \times 10^{-10}$, $\hat{E}_{\max} = 6.3574 \times 10^{-10} \implies p = 4.01$
实测平均收敛阶达 **$4.0191$**，自适应求解器 `simulate_step_adaptive` 通过局部误差 PI 控制律动态调节 $h_k$，在非刚性光滑系统（$S \le 100$）上控制局部截断误差。

### 3. 缩放平方矩阵指数、开环复射线稳定半径与采样闭环谱半径分析
1. **缩放平方矩阵指数 (`matrix_exponential` & `simulate_step_zoh_expm`)**：
   对低阶状态空间矩阵（$n \le 6$），选取缩放幂次 $2^s$ 使 $\|\mathbf{A} dt / 2^s\|_\infty \le 0.5$，采用 12 阶泰勒级数计算 $\exp(\mathbf{A} dt / 2^s)$ 后执行 $s$ 次平方还原，消除线性 ZOH 离散递推的常微分方程步进截断误差。
2. **开环 LTI 二维复射线方向稳定上限**：
   仅当全体极点严格满足 Hurwitz 条件 $\text{Re}(\lambda_i) < 0$ 时开环系统存在正稳定步长（若存在 $\text{Re}(\lambda_i) \ge 0$，则 $h_{\max,\text{open}} = 0$）。对每个左半平面极点 $\lambda_i = |\lambda_i| e^{j\theta_i}$，沿射线 $z = r e^{j\theta_i}$ 二分求解 $|R(r e^{j\theta_i})| = 1$ 得到方向稳定半径 $r_{\text{crit}}(\theta_i)$，从而开环稳定步长上限为：
   $$h_{\max,\text{open}} = \min_{1 \le i \le n} \frac{r_{\text{crit}}(\arg \lambda_i)}{|\lambda_i|}$$
3. **采样数据数字 PID 闭环特征多项式谱半径分析**：
   在离散采样 PID + 零阶保持器（ZOH）线性闭环中，对增广状态 $\mathbf{z}_k = [\mathbf{x}_k^T, I_{k-1}, e_{k-1}]^T \in \mathbb{R}^{n+2}$ 构造单步转移矩阵 $\mathbf{\Phi}_{\text{CL}}(h)$，通过 Faddeev-LeVerrier 递推提取特征多项式 $\det(z\mathbf{I} - \mathbf{\Phi}_{\text{CL}}(h))$ 并求解特征根模长最大值 $\rho(\mathbf{\Phi}_{\text{CL}}(h)) = \max_i |z_i|$。该分析适用于低阶良态系统（$n \le 6$，特征根偏离单位圆 $||z_i| - 1| > 10^{-6}$），可有效识别大采样周期 ZOH 延迟失稳与未滤波高微分增益失稳。
