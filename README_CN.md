# PromPredict Toolkit

> **一键可视化 + 质控 + 文档。让 PromPredict 从命令行工具变研究成果。**  
> 基于 DNA 双链热力学稳定性的启动子预测一站式工具包。

[![Python](https://img.shields.io/badge/python-3.7%2B-blue)](https://www.python.org/)
[![License](https://img.shields.io/badge/license-MIT-green)](LICENSE)
[![DOI](https://img.shields.io/badge/DOI-10.1039%2Fb906535k-red)](https://doi.org/10.1039/B906535K)
[![Platform](https://img.shields.io/badge/platform-linux%20%7C%20wsl-lightgrey)]()

**GitHub:** [github.com/lijiangyong314/PromPredict-Toolkit](https://github.com/lijiangyong314/PromPredict-Toolkit)

---

## 工具包简介

**PromPredict Toolkit** 围绕 IISc（印度科学理工学院）开发的 PromPredict 启动子预测工具，提供了一套完整的科研辅助解决方案。

做启动子预测时，你通常会遇到这些问题：

- ❌ PromPredict 输出了三个文本文件，全是数字，看不懂
- ❌ 几十条序列跑出来上百个预测，哪些可信？哪些是噪音？
- ❌ 审稿人要图——你得自己手搓 matplotlib
- ❌ 基因组模式和多序列模式输出格式不一样，脚本报错找不到问题
- ❌ 想把 PromPredict 和 iPro-MP、promoterCalculator 联合使用，但不知道怎么串起来

这个工具包解决的就是这些事：

| 你需要的 | 工具包提供 |
|----------|-----------|
| **看懂结果** | 自由能剖面图 + 置信度柱状图，一条序列一张 PNG，Dmax 自动标注 |
| **筛选高可信预测** | `awk` 一行命令按 Dmax 阈值过滤，绿/橙/红三色可视化 |
| **质量控制** | 覆盖率检查、Dmax 分布统计、异常值检测脚本 |
| **出版级图表** | `plot_prompredict.py`，200 DPI，可直接投稿 |
| **格式兼容** | 多序列模式 + 全基因组模式双支持，列索引差异文档化 |
| **工具链整合** | 与 iPro-MP、promoterCalculator、ViennaRNA 的联合使用指南 |

**一句话：** 下载 PromPredict 二进制 → 跑数据 → 用本工具包可视化 + 质控 → 得到能直接放进论文的结果。

---

## 功能特性

### 核心：`plot_prompredict.py` 可视化引擎

```bash
python3 plot_prompredict.py sample_seq
```

一条命令生成：

```
sample_seq_profiles/
├── sample_seq_confidence_summary.png   ← 所有预测按 Dmax 排序的汇总柱状图
├── PM0-4641yaiS...png                  ← 每条序列独立一张自由能剖面图
├── PM00066ompF...png                   ← 蓝线=ΔG曲线，红带=预测启动子
└── ...
```

每张剖面图包含：ΔG 折线、全序列均值线、预测启动子色带（绿/橙/红三色）、Dmax 数值标注、Dmax-ΔG 一致性警告。

### 内置质控管线

```bash
# 覆盖率 | Dmax 分布 | 异常值 | 按置信度分级统计
# 全部内置在 README 中的复制即用命令
```

### 双模式支持

| 模式 | 典型场景 | 可视化 |
|------|---------|--------|
| 多序列模式 | 几十条基因上游区逐一预测 | ✅ `plot_prompredict.py` 直接支持 |
| 全基因组模式 | 扫描完整染色体 | ⚠️ 命令行筛选 + 计划中的浏览器轨迹图 |

### 完整的操作手册（这份 README）

涵盖安装、输入/输出格式、QC、参数调优决策树、与 AI 工具联合使用、适用范围、FAQ。

---

---

## 示例输出

光说不练假把式。`sample/` 目录里放的是用 IISc 官方样本数据集跑出来的真实结果——输入、输出、可视化全在里面。

### `sample/data/` 里面有什么

| 文件 | 内容 |
|------|------|
| `sample_seq_stb.txt` | 自由能（ΔG）剖面：10 个已知 *E. coli* 启动子区的 15-bp 滑动窗口数据 |
| `sample_seq_PPde.txt` | 预测启动子列表：14 个预测，Dmax 从 1.02（噪音）到 4.54（强信号） |
| `sample_seq_GCstat.txt` | 输入序列的 GC 含量分布 |

### `sample/profiles/` 里面有什么

由 `plot_prompredict.py` 生成——**你跑自己数据后拿到的就是这种图**：

| 文件 | 展示内容 |
|------|---------|
| `sample_seq_confidence_summary.png` | 14 个预测按 Dmax 排序。绿 = 高（≥3.0），橙 = 中（2.0–3.0），红 = 低（<2.0）。一眼看出：4 个强候选、6 个中等、4 个弱信号 |
| `PM0-*.png` | 每条序列一张自由能剖面图。蓝线 = ΔG，红色带 = 预测启动子，Dmax 数值直接标在图上。经典的大肠杆菌 ompF 启动子被正确识别（Dmax=2.34） |

### 自己试试

```bash
# 原始输入在 sample/data/ 里，跑一遍：
printf "sample/data/sample_seq.txt\n100\ndefault\n" | ./prompredict

# 可视化：
python3 plot_prompredict.py sample_seq

# 你的输出应该和 sample/profiles/ 一致
```

> 💡 IISc 样本数据集覆盖 10 个实验验证的 *E. coli* σ⁷⁰ 启动子（yaiS、ompF、pgaA、csgD、csgB、rpoB、gerD、skf、glpT、srfAA），PromPredict 全部正确识别。

---

## 快速开始

```bash
# 1. 克隆仓库
git clone https://github.com/lijiangyong314/PromPredict-Toolkit.git
cd PromPredict-Toolkit

# 2. 安装 Python 依赖（可视化用）
pip install matplotlib numpy

# 3. 下载 PromPredict 二进制
wget --no-check-certificate \
  "https://dna.mbu.iisc.ac.in/prompredict/exe/linux/PromPredict_mulseq" \
  -O prompredict && chmod +x prompredict

# 4. 跑数据
printf "test.fasta\n100\ndefault\n" | ./prompredict

# 5. 可视化
python3 plot_prompredict.py test
```

---

## 安装

### 前置条件

- **Python** 3.7+（仅可视化需要）
- **matplotlib + numpy**（`pip install matplotlib numpy`）
- **操作系统**：Linux x86-64 或 WSL2
- PromPredict 二进制约 2.4 MB

### 获取 PromPredict 二进制

```bash
# 多序列模式（最常用）
wget --no-check-certificate \
  "https://dna.mbu.iisc.ac.in/prompredict/exe/linux/PromPredict_mulseq" \
  -O prompredict && chmod +x prompredict

# 全基因组模式（>10 MB 基因组用）
wget --no-check-certificate \
  "https://dna.mbu.iisc.ac.in/prompredict/exe/linux/PromPredict_genome_V2" \
  -O prompredict_genome && chmod +x prompredict
```

> IISc 服务器 SSL 证书已过期，**必须**加 `--no-check-certificate`。

### Conda 隔离环境（可选）

```bash
conda create -n prompredict python=3.9 matplotlib numpy -y
conda activate prompredict
```

---

## 底层工具：PromPredict

本工具包依赖于 IISc 开发的 **PromPredict** 二进制文件。以下是其核心信息。

### 原理

PromPredict 基于一个物理直觉：**启动子区域的 DNA 必须是局部不稳定的，RNA 聚合酶才能将双链解开并启动转录。** 它使用 Breslauer 二核苷酸自由能参数，对输入序列上每一个 15-bp 滑动窗口计算平均自由能（ΔG）。当某个区域的 ΔG 值显著高于两侧基线时，即被标记为候选启动子。

| 属性 | 值 |
|------|-----|
| **来源** | 印度科学理工学院 (IISc Bangalore), Molecular Biophysics Unit |
| **论文** | *Molecular BioSystems* (2009); *Scientific Reports* (2018) |
| **方法** | 热力学 (ΔG)，非机器学习 |
| **依赖** | 无（独立二进制） |
| **大小** | 2.4 MB |
| **训练/阈值** | *E. coli*、*B. subtilis*、*M. tuberculosis* 实验验证 TSS |

### 与其他工具对比

| 方法 | 原理 | 模型大小 | 需 GPU | 可解释性 |
|------|------|---------|--------|---------|
| **PromPredict** | 热力学 (ΔG) | 2.4 MB 二进制 | 否 | ✅ 每个预测可推理 |
| iPro-MP | BERT Transformer | 36.5 GB | 是 | ❌ |
| Promotech | 随机森林 / RNN | 自定义 | 可选 | 部分 |

### 适用范围

| 场景 | 适用度 |
|------|--------|
| 革兰氏阴性菌（变形菌门），GC 35–55% | ✅ 最佳 |
| 其他细菌、低等真核（酵母、线虫、拟南芥） | ✅ 可用 |
| 高 GC 细菌（*Streptomyces* >72%）或古菌 | ⚠️ 需手动调参 |
| 哺乳动物基因组 | ❌ 不推荐 |
| 预测启动子**强度**（而非位置） | ❌ 需用 promoterCalculator |

> ⚠️ **不确定性声明**：Dmax 置信度分级（≥3.5 = 极高，3.0–3.5 = 高，2.0–3.0 = 中等，<2.0 = 低）是我们在 IISc 样本数据集和文献验证的 *E. coli* 启动子上测试得出的经验性指南，并非 PromPredict 原论文中的硬性标准。在不同物种或极端 GC 条件下，最优截断值可能偏移。建议在广泛使用前，先用目标物种的已知启动子进行基准测试。

---

## 使用方法

### 交互模式

直接运行二进制，依次回答三个参数：

```
Enter the Input File Name:        your_file.fasta
Enter the E1 region window size:  100
Whole genome GC content:          default
```

### 非交互模式（推荐）

```bash
printf "input.fasta\n100\ndefault\n" | ./prompredict
```

`printf` 通过管道自动填入三行答案。

| 参数 | 默认值 | 含义 |
|------|--------|------|
| 输入文件 | *(必填)* | 多 FASTA 文件 |
| E1 窗口 | `100` | 上游比较窗口 (nt) |
| GC 含量 | `default` | `default`=自动检测；也可填数字如 `50.5` |

### 全基因组模式

```bash
printf "genome.fna\n100\ndefault\n" | ./prompredict_genome
```

> ⚠️ 全基因组模式输出格式与多序列模式不同（11 列 vs. 8 列），见[输出格式详解](#输出格式详解)。

---

## 输出格式详解

### 多序列模式

运行后生成三个文件（以 `input.fasta` 为例）：

#### `input_PPde.txt` — 预测启动子列表

每个启动子两行：

```
ID    <序列ID>    <GC%>
><起始>..<终止>    <长度>    <序列>    <lsp>    <lspe>    <Dmax位置>    <Dmax>    <Dave>
```

| 列 | 含义 | 判断标准 |
|----|------|---------|
| `>起始..终止` | 启动子坐标（**1-based**） | — |
| `lsp` | 局部稳定性峰值，≈ −10 区 | — |
| **`Dmax`** | **核心：启动子与侧翼的 ΔG 差值** | ≥3.5 极高 / 3.0-3.5 高 / 2.0-3.0 中 / <2.0 低 |
| `Dave` | 平均稳定性差 | >1.5 较好 |

> 💡 PromPredict 用 1-based 坐标。Python 提取序列时：`start_py = start - 1`。

#### `input_stb.txt` — 自由能剖面

```
ID    <序列ID>    ΔG₁    ΔG₂    ΔG₃    ...
```

越负（−22）= 越稳定 → 不是启动子；越正（−12）= 不稳定 → 候选启动子。行数 = 序列长度 − 14。

#### `input_GCstat.txt` — GC 含量分布

内部使用，通常无需查看。

### 全基因组模式（格式差异）

`_PPde.txt` 为扁平 11 列格式（无 ID/>标题行）：

```
# 窗口起始  %GC  pstart  pend  长度  序列  lsp  lspe  Dmax位置  Dmax  Dave
0    15.90   500   607   108   tcgcaaat...   549   -11.42   555   1.85   1.52
```

| 关键差异 | 多序列模式 | 全基因组模式 |
|---------|-----------|-------------|
| Dmax 列索引 | `$7` | **`$10`** |
| Dave 列索引 | `$8` | **`$11`** |
| pstart/pend | `parts[0]` 解析 | `$3` / `$4` |

---

## 质量控制

### 多序列模式

```bash
# 覆盖率
grep -c "^>" input.fasta        # 输入序列数
grep -c "^ID" input_PPde.txt    # 检出 ≥1 个启动子的序列

# Dmax 分布
awk -F'\t' '/^>/ {print $7}' input_PPde.txt | sort -nr | head -20

# 按置信度分级
high=$(awk -F'\t' '/^>/ && $7>=3.0' input_PPde.txt | wc -l)
med=$(awk  -F'\t' '/^>/ && $7>=2.0 && $7<3.0' input_PPde.txt | wc -l)
low=$(awk  -F'\t' '/^>/ && $7<2.0' input_PPde.txt | wc -l)
echo "高: $high | 中: $med | 低: $low"
```

若全部 Dmax < 2.0：① 检查序列长度 ≥ 200 bp ② 高 GC 基因组指定实际 GC% ③ 改用 E1=50。

### 全基因组模式

```bash
# 注意列索引不同
awk -F'\t' 'NR>2 && $10 >= 3.0' genome_PPde.txt | wc -l         # 高
awk -F'\t' 'NR>2 && $10 >= 2.0 && $10 < 3.0' genome_PPde.txt | wc -l  # 中
awk -F'\t' 'NR>2 && $10 < 2.0' genome_PPde.txt | wc -l          # 低

# 导出高置信坐标（BED 格式）
awk -F'\t' 'NR>2 && $10 >= 3.0 {print $1"\t"$4"\t"$5"\t"$9"\t"$10}' genome_PPde.txt > high_conf.bed
```

---

## 可视化：plot_prompredict.py

```bash
python3 plot_prompredict.py <prefix>
```

输出目录 `<prefix>_profiles/`：

| 文件 | 内容 |
|------|------|
| `<prefix>_confidence_summary.png` | 所有预测按 Dmax 降序的横向柱状图 |
| `<序列ID>.png` | 每条序列一张 16×5 宽幅剖面图 |

每张剖面图包含：

- 🔵 蓝色折线 — 15-bp 窗口 ΔG 曲线
- 🔴 红色虚线 — 序列平均 ΔG
- 🟢/🟠/🔴 色带 — 预测启动子（绿=高/橙=中/红=低置信度）
- 数值标注 — Dmax 值标注在色带上方
- ⚠ 警告 — Dmax 与 ΔG 不一致时标紫色警告

> ⚠️ 目前仅支持多序列模式。全基因组可视化（浏览器风格轨迹图）正在开发中。

---

## 参数调优（决策树）

| 现象 | 原因 | 方案 |
|------|------|------|
| 假阳性多 | 阈值过低 | Dmax 截断值提高到 3.5 |
| 假阴性多 | 100nt 窗口掩盖短信号 | E1=50 重跑 |
| 高 GC（>60%） | ΔG 基线偏低 | 指定实际 GC%，放宽到 ≥2.0 |
| 低 GC（<30%） | ΔG 基线偏高 | 指定实际 GC% |

---

## 与其他工具联合使用

PromPredict 做**初筛**，结合以下工具构建完整注释管线：

| 步骤 | 工具 | 互补点 |
|------|------|--------|
| **初筛** | PromPredict | 唯一物理原理方法，解释"为什么是启动子" |
| **验证** | iPro-MP / DeePromoter | 识别 σ 因子特异性 motif |
| **强度** | promoterCalculator | 量化 −35/−10 匹配度 |
| **终止** | cryptkeeper / transtermHP | 启动子常与终止子共定位 |
| **结构** | ViennaRNA | RNA 二级结构影响转录效率 |

---

## 常见问题

<details>
<summary><b>_stb.txt 行数为什么比序列长度少？</b></summary>
15-bp 窗口每端损失 7 bp，行数 = 序列长度 − 14。
</details>

<details>
<summary><b>预测的启动子在 CDS 内部？</b></summary>
不一定错误。可能存在内部启动子、终止子附近解链区。结合 −10/−35 盒 motif 验证。
</details>

<details>
<summary><b>每条序列都检出启动子？</b></summary>
序列 >500 bp 时首轮较宽松。按 Dmax ≥ 2.0 筛选即可。
</details>

<details>
<summary><b>macOS 能用吗？</b></summary>
二进制为 Linux x86-64。用 Docker 或 IISc 在线版。
</details>

<details>
<summary><b>不确定：Dmax 阈值是通用的吗？</b></summary>
推荐的分级标准（≥3.5, 3.0–3.5, 2.0–3.0, <2.0）是在 IISc 样本数据和已知 <i>E. coli</i> 启动子上测试得到的经验性指南。不同物种或极端 GC 条件下最优阈值可能偏移。建议先用目标物种的几个已知启动子做基准测试，再全基因组应用。
</details>

---

## 引用

> Rangannan, V., & Bansal, M. (2009).  
> **Relative stability of DNA as a generic criterion for promoter prediction.**  
> *Molecular BioSystems*, 5(12), 1758–1769. DOI: [10.1039/B906535K](https://doi.org/10.1039/B906535K)

真核生物应用另引：

> Rangannan, V., & Bansal, M. (2018).  
> **Identification of putative promoters in 48 eukaryotic genomes.**  
> *Scientific Reports*, 8, 4520. DOI: [10.1038/s41598-018-22129-8](https://doi.org/10.1038/s41598-018-22129-8)

---

## 贡献

欢迎 Issues 和 PR：[github.com/lijiangyong314/PromPredict-Toolkit](https://github.com/lijiangyong314/PromPredict-Toolkit)

## 许可证

本工具包（Python 脚本、文档、可视化代码）以 **MIT 许可证**发布。详见 [LICENSE](LICENSE)。

PromPredict 二进制 © IISc Molecular Biophysics Unit，分发条款见[官网](https://dna.mbu.iisc.ac.in/prompredict/download.html)。

> ⚠️ **不确定性声明**：PromPredict 二进制文件的具体许可/分发条款以 IISc 官网为准。以上版权归属基于网站页脚标注，请查阅官网获取确切的法律文本。

## 致谢

- Vetriselvi Rangannan & Manju Bansal（IISc）
- Breslauer et al. (1986 *PNAS*) 二核苷酸参数
- *E. coli* K-12、*B. subtilis*、*M. tuberculosis* 实验 TSS 数据集
