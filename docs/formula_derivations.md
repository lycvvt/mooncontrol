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

在 `src/sim/sim_wbtest.mbt` 中，数值仿真 RK4 的计算结果与上述解析公式的最大绝对误差 $< 10^{-14}$。

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

RK4 具备四阶局部截断误差 $O(dt^5)$ 与四阶整体截断误差 $O(dt^4)$，在步长 $dt = 0.0005$ 时仿真精度达 $10^{-14}$ 量级。
