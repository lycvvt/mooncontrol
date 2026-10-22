# MoonControl 🌙🕹️

> **MoonBit 上的控制系统设计与闭环仿真工具链**  
> Control System Design and Simulation Toolchain in Pure MoonBit (Zero FFI)

[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](LICENSE)
[![MoonBit](https://img.shields.io/badge/MoonBit-0.1.20260920-purple.svg)](https://www.moonbitlang.com/)
[![Tests](https://img.shields.io/badge/Tests-29%20passed-brightgreen.svg)]()
[![Hackathon](https://img.shields.io/badge/MoonBit%20Hackathon-2026%20October-green.svg)](https://moonbitlang.github.io/Hackathon2026/)

---

## 📌 生态定位与背景

截至 2026 年 10 月，在 Mooncakes 官方全量中心（2,800+ 模块）与 GitHub 生态中检索控制理论相关主题：
- 检索 `PID / control system / transfer function / 状态空间 / Bode / root locus / Lyapunov`
- Mooncakes 仅有状态估计类库（如 Kalman 滤波）与形式化验证库，**没有任何面向控制器设计、经典控制理论与闭环仿真的工具链**。
- GitHub 检索 `moonbit control system` 与 `moonbit pid` 结果均为 0。

**MoonControl** 旨在填补这一生态空白，为 MoonBit 开发者与系统工程师提供一套纯 MoonBit 实现、**完全零 FFI** 的高精度控制系统分析、设计与时域仿真基础库。

---

## ✨ 核心特性

- **纯原生零 FFI (Zero FFI)**：复数四则运算、多项式代数、Durand-Kerner 根求解与高阶微分方程积分器均使用纯 MoonBit 原生代码编写，不依赖任何 C/JS 外部库。
- **高精度解析解基准**：严格以一阶、二阶系统（欠阻尼/临界阻尼/过阻尼）解析阶跃响应为数学基准，数值仿真误差收敛于 `< 1e-14`（远超规范要求的 `< 1e-9`）。
- **经典控制系统工具箱**：
  - **传递函数 (Transfer Function)**：多项式分子分母建模，支持串联 (Series)、并联 (Parallel)、反馈连接 (Feedback)，自动提取极点与零点。
  - **时域分析与性能指标**：超调量 $M_p$、峰值时间 $t_p$、上升时间 $t_r$、调节时间 $t_s$（2% 与 5% 判据）、终值定理稳态误差 $e_{ss}$。
  - **稳定性分析**：Routh-Hurwitz 代数判据（支持首列元素为 0 微扰与全零行辅助多项式求导退化处理），与极点实部联合判定。
  - **PID 控制器**：连续形式与离散位置式 PID，内置积分限幅（Anti-windup）与一阶低通微分滤波。
  - **时域仿真引擎**：经典四阶 Runge-Kutta (RK4) 高精度数值积分器。
  - **可视化图表**：自研轻量级 SVG 矢量图与终端 ASCII 字符曲线生成，零外部绘图依赖。

---

## 🚀 快速开始 (Quick Start)

### 1. 运行完整测试套件
```bash
moon test
```
全部 29 个测试用例秒级通过，覆盖 8 组控制理论基准验证。

### 2. 运行主程序 CLI 交互演示
```bash
moon run cmd/main
```
输出包含：传递函数建模、极点提取、一阶解析解基准核对、RK4 误差断言、终端 ASCII 阶跃曲线、二阶 PID 闭环仿真与 SVG 生成。

### 3. 运行精选示例
```bash
# 示例 1: 一阶系统阶跃响应与解析解对比
moon run examples/first_order

# 示例 2: 二阶被控对象闭环离散 PID 整定与时域响应
moon run examples/second_order_pid

# 示例 3: Routh-Hurwitz 阵列与零极点分布稳定性分析
moon run examples/stability_analysis
```

---

## 📊 可视化效果展示

### 1. 终端 ASCII 阶跃响应曲线 (Terminal ASCII Curve)
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

### 2. 原生矢量图输出 (SVG Vector Graphic)
通过 `@plot.plot_svg` 自动生成包含坐标网格、时间轴刻度、设定值参考线（Setpoint）与抗锯齿响应折线的独立矢量图文件。

---

## 🧪 8 组核心验证测试矩阵 (Definition of Done)

根据赛事技术规范，本项目实现了严密对齐控制理论公式的 8 组基准测试：

| # | 测试组 | 验证内容与数学基准 | 误差指标 | 状态 |
|---|---|---|---|---|
| **1** | **复数代数** | 复数加减乘除、模长、幅角、极坐标与欧拉公式 $e^{i\pi} = -1$ | 手算精确对齐 |  PASSED |
| **2** | **多项式算法** | Horner 求值、多项式乘法卷积、Durand-Kerner 复根求解 | 代数精确对齐 |  PASSED |
| **3** | **一阶阶跃响应** | $G(s) = \frac{K}{Ts+1}$ 在 $t=T, 2T, 3T$ 与解析解 $K(1-e^{-t/T})$ 比对 | **$7.7 \times 10^{-15} < 10^{-9}$** |  PASSED |
| **4** | **二阶阶跃响应** | 欠阻尼 ($0<\zeta<1$)、临界阻尼 ($\zeta=1$)、过阻尼 ($\zeta>1$) 解析解 | **$< 10^{-14} < 10^{-9}$** |  PASSED |
| **5** | **动态性能指标** | 超调量 $M_p$、峰值时间 $t_p$、调节时间 $t_s$ 公式验算 | **$< 10^{-7} < 10^{-6}$** |  PASSED |
| **6** | **稳态误差分析** | 0 型/1 型系统在阶跃与斜坡输入下，对比拉氏变换终值定理 | 绝对一致 |  PASSED |
| **7** | **稳定性交叉验证** | Routh-Hurwitz 阵列判定与极点实部分布 $\text{Re}(p) < 0$ 严格对齐 | 符号一致性 100% |  PASSED |
| **8** | **闭环时域仿真** | 受控对象闭环仿真轨迹与标准二阶解析解严格对齐 | **$1.6 \times 10^{-14} < 10^{-6}$** |  PASSED |

---

## 🏗️ 架构与模块划分

```
mooncontrol/
├── moon.mod                    # MoonBit 模块清单
├── LICENSE                     # Apache-2.0 许可证
├── README.md                   # 项目介绍与快速开始
├── docs/
│   └── formula_derivations.md  # 详尽控制理论公式推导与算法笔记
├── src/
│   ├── num/                    # 数值基础：复数、多项式 (Horner / Durand-Kerner)、矩阵
│   ├── tf/                     # 传递函数：建模、串联、并联、反馈闭环、零极点与直流增益
│   ├── analysis/               # 解析解引擎、动态性能指标、稳态误差、Routh-Hurwitz 判据
│   ├── pid/                    # 连续/离散 PID 控制器（带 Anti-windup 积分限幅与滤波）
│   ├── sim/                    # 可控标准型状态空间实现与经典四阶 Runge-Kutta (RK4) 仿真
│   └── plot/                   # 终端 ASCII 字符图与 SVG 矢量图渲染器
├── cmd/
│   └── main/                   # 交互式 CLI 演示工具
└── examples/
    ├── first_order/            # 示例 1: 一阶系统解析阶跃响应
    ├── second_order_pid/       # 示例 2: 二阶系统与离散 PID 闭环追踪
    └── stability_analysis/     # 示例 3: Routh-Hurwitz 与极点稳定性判定
```

---

## 🚫 非目标 (Non-Goals)

为保持工程的轻量与核心质量，本项目明确**不包含**以下特性：
- 非线性控制、鲁棒控制（$H_\infty$）、模型预测控制 (MPC)
- 真实硬件在环 (HIL) 通信驱动
- 重量级图形界面 (GUI) 框架依赖
- 面向 PLC 的自动代码生成

---

## 📚 理论依据与参考来源

本项目相关数学模型与算法严格对照经典控制理论标准教材推导验证：
1. Katsuhiko Ogata, *Modern Control Engineering* (5th Edition), Pearson.
2. Richard C. Dorf & Robert H. Bishop, *Modern Control Systems* (13th Edition), Pearson.
3. [python-control](https://github.com/python-control/python-control) - MATLAB/Simulink 接口设计与测试用例对照.
4. 详细数学推导请查阅 [公式推导笔记](docs/formula_derivations.md).

---

## 📄 许可证

本项目采用 [Apache-2.0 License](LICENSE) 开源。
