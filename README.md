# MoonControl

> **基于 MoonBit 实现的控制系统建模、分析与高精度仿真工具链 (Zero FFI)**  
> A pure MoonBit toolchain for linear control systems design, stability analysis, and time-domain simulation.

[![CI](https://github.com/lycvvt/mooncontrol/actions/workflows/ci.yml/badge.svg)](https://github.com/lycvvt/mooncontrol/actions)
[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](LICENSE)
[![Language](https://img.shields.io/badge/Language-MoonBit-purple.svg)](https://www.moonbitlang.com/)
[![Tests](https://img.shields.io/badge/Tests-29%2F29%20passed-success.svg)]()

---

## 项目概述

MoonControl 是首个为 MoonBit 生态构建的原生控制理论基础库与时域仿真系统，提供连续时间线性时不变（LTI）系统建模、代数稳定性判据、动态性能指标计算、连续/离散 PID 控制以及经典四阶 Runge-Kutta（RK4）数值积分。

### 生态定位
截至 2026 年 10 月，在 Mooncakes 官方中心（2,800+ 模块）与开源社区中：
- 检索 `control system`、`transfer function`、`PID`、`Routh-Hurwitz`、`state space`，除状态估计（Kalman 滤波）外，控制系统建模与闭环仿真领域尚无成熟工具链。
- GitHub 检索 `moonbit control system` 与 `moonbit pid` 结果均为 0。

MoonControl 填补了该基础工程层级，面向嵌入式控制、系统工程与高校自动化教学提供高性能、零 FFI 依赖的确定性计算基础设施。

### 核心指标
- **零 FFI (Zero Foreign Function Interface)**：复数代数、多项式根求解、矩阵运算、积分引擎全部使用纯 MoonBit（`.mbt`）原生实现。
- **解析解数值对齐**：以控制理论一阶、二阶系统阶跃响应解析闭式解为基准，数值仿真全轨迹误差控制在 $10^{-14}$ 量级（严格优于规范设定的 $10^{-9}$ 要求）。
- **确定性与无异常崩溃**：核心分析与仿真函数均为纯函数设计，非法输入统一返回 `Result[T, String]`，杜绝运行时非预期 panic。

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

## 理论验证与测试基准

本项目严格建立 8 组对齐经典控制理论基准公式的自动化测试套件：

| 编号 | 测试模块 | 验证基准与技术点 | 误差指标 | 判定 |
| :---: | :--- | :--- | :---: | :---: |
| **1** | **复数运算** | 加减乘除、共轭、模长、欧拉公式 $e^{i\pi} = -1$ | 手算精确对齐 | PASS |
| **2** | **多项式代数** | Horner 求值、多项式乘法卷积、Durand-Kerner 根求解 | 代数精确对齐 | PASS |
| **3** | **一阶阶跃响应** | $G(s) = \frac{K}{Ts+1}$ 在 $t=T, 2T, 3T$ 与解析式 $K(1-e^{-t/T})$ 比对 | **$7.7 \times 10^{-15} < 10^{-9}$** | PASS |
| **4** | **二阶阶跃响应** | 欠阻尼 ($0<\zeta<1$)、临界阻尼 ($\zeta=1$)、过阻尼 ($\zeta>1$) 解析解 | **$< 10^{-14} < 10^{-9}$** | PASS |
| **5** | **动态性能指标** | 欠阻尼系统超调量 $M_p$、峰值时间 $t_p$、调节时间 $t_s$ 公式验算 | **$< 10^{-7} < 10^{-6}$** | PASS |
| **6** | **稳态误差分析** | 0 型/1 型系统在阶跃与斜坡输入下对照拉氏变换终值定理 | 理论完全一致 | PASS |
| **7** | **稳定性交叉检验** | Routh-Hurwitz 阵列符号变化与极点实部分布 $\text{Re}(p) < 0$ 比对 | 判定一致率 100% | PASS |
| **8** | **闭环时域仿真** | 受控对象闭环 RK4 仿真轨迹与标准二阶解析解全时域比对 | **$1.6 \times 10^{-14} < 10^{-6}$** | PASS |

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
