# MoonControl 🌙🕹️

> **MoonBit 上的控制系统设计与仿真工具链**  
> Control System Design and Simulation Toolchain in Pure MoonBit (Zero FFI)

[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](LICENSE)
[![MoonBit](https://img.shields.io/badge/MoonBit-0.1.20260920-purple.svg)](https://www.moonbitlang.com/)
[![Hackathon](https://img.shields.io/badge/MoonBit%20Hackathon-2026%20October-green.svg)](https://moonbitlang.github.io/Hackathon2026/)

---

## 📌 生态定位与背景

截至 2026 年 10 月，在 Mooncakes 官方全量中心（2,800+ 模块）与 GitHub 生态中检索控制理论相关主题：
- `PID / control system / transfer function / 状态空间 / Bode / root locus / Lyapunov`
- Mooncakes 仅有状态估计类库（如 Kalman 滤波）与形式化验证库，**没有任何面向控制器设计、经典控制理论与闭环仿真的工具链**。
- GitHub 检索 `moonbit control system` 与 `moonbit pid` 结果均为 0。

**MoonControl** 旨在填补这一生态空白，为 MoonBit 开发者与嵌入式/系统工程师提供一套纯 MoonBit 实现、**完全零 FFI** 的高精度控制系统分析、设计与时域仿真基础库。

---

## ✨ 核心特性

- **纯原生零 FFI (Zero FFI)**：复数四则运算、多项式代数、根求解与高阶微分方程积分器均使用纯 MoonBit 原生代码编写，不依赖任何 C/JS 外部库。
- **高精度解析解基准**：严格以一阶、二阶系统（欠阻尼/临界阻尼/过阻尼）解析阶跃响应为数学基准，数值仿真误差收敛于 `< 1e-9`。
- **经典控制系统工具箱**：
  - **传递函数 (Transfer Function)**：多项式分子分母建模，支持串联 (Series)、并联 (Parallel)、反馈连接 (Feedback)，自动提取极点与零点。
  - **时域分析与性能指标**：超调量 $M_p$、峰值时间 $t_p$、上升时间 $t_r$、调节时间 $t_s$（2% 判据）、终值定理稳态误差 $e_{ss}$。
  - **稳定性分析**：Routh-Hurwitz 代数判据（支持首列元素为 0 扰动与全零行辅助多项式退化处理），与极点实部联合判定。
  - **PID 控制器**：连续形式与离散位置式 PID，内置积分限幅（Anti-windup）与一阶低通微分滤波。
  - **时域仿真引擎**：经典四阶 Runge-Kutta (RK4) 高精度数值积分器。
  - **可视化图表**：自研轻量级 SVG 矢量图与终端 ASCII 字符曲线生成，零外部绘图依赖。

---

## 🏗️ 架构与模块划分

```
mooncontrol/
├── moon.mod               # 模块清单
├── src/
│   ├── num/               # 基础数学库：复数、多项式（Horner/Durand-Kerner）、矩阵
│   ├── tf/                # 传递函数：建模、串并联与反馈联结、零极点
│   ├── analysis/          # 控制理论解析解、性能指标计算、Routh-Hurwitz 稳定性判据
│   ├── sim/               # 闭环时域仿真引擎（RK4 积分器）
│   ├── pid/               # 连续/离散 PID 控制器与抗饱和算法
│   └── plot/              # SVG 矢量折线图与 ASCII 终端绘图
└── cmd/
    └── main/              # 命令行交互工具（CLI）
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

---

## 📄 许可证

本项目采用 [Apache-2.0 License](LICENSE) 开源。
