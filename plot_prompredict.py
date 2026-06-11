#!/usr/bin/env python3
"""
PromPredict 自由能剖面可视化 v2.0
用法：python3 plot_prompredict.py sample_seq
输出：sample_seq_profiles/ 文件夹，每序列一张图 + 汇总柱状图

改动日志：
  v2.0 - dataclass 重构 + 防御性IO + 出版级绘图 + 生物学校验
  v1.0 - 初始版本
"""

import sys, os, argparse, logging
from dataclasses import dataclass, field
from collections import defaultdict
from typing import List, Dict, Tuple

import matplotlib.pyplot as plt
import numpy as np

# ====== 0. 科学与视觉常量 ======

THRESHOLDS = {
    "min_windows": 10,       # <10 窗口的序列剔除（<150nt 无意义）
    "dmax_high": 3.0,        # 高置信阈值
    "dmax_medium": 2.0,      # 中等置信阈值
}

COLORS = {
    "stability_baseline":  "#4A90D9",   # ΔG 折线
    "stability_unstable":  "#FF6B6B",   # 不稳定区 = 启动子信号
    "confidence_high":     "#27AE60",   # Dmax >= 3.0
    "confidence_medium":   "#F39C12",   # Dmax 2.0-3.0
    "confidence_low":      "#E74C3C",   # Dmax < 2.0
    "border_promoter":     "#CC0000",   # 启动子色带边框
    "mismatch_warning":    "#9B59B6",   # Dmax-ΔG 不一致警告
}

logging.basicConfig(
    level=logging.INFO,
    format="%(asctime)s [%(levelname)s] %(message)s",
    datefmt="%H:%M:%S",
)

# ====== 1. 数据结构 ======

@dataclass
class StabilityProfile:
    """一条序列的自由能剖面"""
    seq_id: str
    positions: List[int] = field(repr=False)
    energies: List[float] = field(repr=False)

    @property
    def length(self) -> int:
        return len(self.positions)

    @property
    def mean_dg(self) -> float:
        return float(np.mean(self.energies))

    @property
    def dg_range(self) -> float:
        return float(max(self.energies) - min(self.energies))

    def region_mean(self, start: int, end: int) -> float:
        """给定区域的 ΔG 均值（0-based 索引）"""
        start = max(0, start)
        end = min(len(self.energies), end)
        if end <= start:
            return self.mean_dg
        return float(np.mean(self.energies[start:end]))


@dataclass
class PromoterRegion:
    """一个预测启动子"""
    start: int
    end: int
    dmax: float
    dave: float

    def overlaps_profile(self, profile: StabilityProfile) -> bool:
        """启动子坐标是否在序列范围内"""
        return 0 < self.start <= profile.length and self.end <= profile.length

    def confidence_label(self) -> str:
        if self.dmax >= THRESHOLDS["dmax_high"]:
            return "High"
        elif self.dmax >= THRESHOLDS["dmax_medium"]:
            return "Medium"
        return "Low"

    def confidence_color(self) -> str:
        if self.dmax >= THRESHOLDS["dmax_high"]:
            return COLORS["confidence_high"]
        elif self.dmax >= THRESHOLDS["dmax_medium"]:
            return COLORS["confidence_medium"]
        return COLORS["confidence_low"]

    def check_dg_agreement(self, profile: StabilityProfile) -> bool:
        """
        校验：Dmax 高 → 该区域 ΔG 应该高于序列均值（更不稳定）
        返回 True = 热力学一致，False = 矛盾
        """
        return profile.region_mean(self.start, self.end) > profile.mean_dg


# ====== 2. IO 与解析 ======

def load_stb(stb_file: str) -> Dict[str, StabilityProfile]:
    """
    解析 PromPredict *_stb.txt 自由能剖面。

    格式：
        # 注释行（跳过）
        ID\tseq_id\tΔG1\tΔG2\t...

    防御：
        - 验证第三列确实是数值（排除表头行误匹配）
        - 拒绝少于 min_windows 个窗口的序列
        - 逐行捕获 ValueError
    """
    profiles = {}
    skipped_empty = 0
    skipped_short = 0
    skipped_parse = 0

    with open(stb_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split("\t")
            if len(parts) < 3:
                skipped_empty += 1
                continue

            # 验证第三列是数字，排除 "stability values for..." 表头行
            if not (parts[0] == "ID" and parts[2].lstrip("-").replace(".", "").isdigit()):
                continue

            seq_id = parts[1]

            try:
                values = [float(v) for v in parts[2:] if v.strip()]
            except ValueError:
                skipped_parse += 1
                continue

            if len(values) < THRESHOLDS["min_windows"]:
                skipped_short += 1
                continue

            positions = list(range(1, len(values) + 1))
            profiles[seq_id] = StabilityProfile(seq_id, positions, values)

    # 汇总统计
    if skipped_short:
        logging.warning(
            f"Skipped {skipped_short} sequence(s) with <{THRESHOLDS['min_windows']} "
            f"windows (too short for meaningful prediction)"
        )
    if skipped_parse:
        logging.warning(f"Skipped {skipped_parse} sequence(s) due to parse errors")
    if not profiles:
        logging.error("No valid profiles loaded. Check input file format.")

    return profiles


def load_ppde(ppde_file: str) -> Dict[str, List[PromoterRegion]]:
    """
    解析 PromPredict *_PPde.txt 预测启动子列表。

    格式（每两行一组）：
        ID\tseq_id\tGC%
        >start..end\tlength\tseq\tlsp\tlspe\tdmax_pos\tdmax\tdave

    防御：
        - current_id 必须在 >行之前被 ID 行赋值
        - 所有数值解析包裹在 try/except 中
    """
    promoters: Dict[str, List[PromoterRegion]] = defaultdict(list)
    current_id = None
    skipped_no_id = 0
    skipped_parse = 0

    with open(ppde_file) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split("\t")

            if parts[0] == "ID":
                current_id = parts[1]

            elif parts[0].startswith(">"):
                if current_id is None:
                    skipped_no_id += 1
                    continue

                try:
                    start, end = map(int, parts[0][1:].split(".."))
                    dmax = float(parts[6])
                    dave = float(parts[7])
                except (ValueError, IndexError):
                    skipped_parse += 1
                    continue

                promoters[current_id].append(
                    PromoterRegion(start, end, dmax, dave)
                )

    if skipped_no_id:
        logging.warning(
            f"Skipped {skipped_no_id} promoter line(s) without preceding ID line. "
            f"Possible file corruption or format mismatch."
        )
    if skipped_parse:
        logging.warning(f"Skipped {skipped_parse} promoter line(s) due to parse errors")

    return promoters


# ====== 3. 视觉常量 ======

def promoter_fill_color(dmax: float) -> str:
    """启动子高亮色带颜色"""
    if dmax >= THRESHOLDS["dmax_high"]:
        return COLORS["confidence_high"]
    elif dmax >= THRESHOLDS["dmax_medium"]:
        return COLORS["confidence_medium"]
    return COLORS["confidence_low"]


# ====== 4. 绘图 ======

def plot_single_profile(
    profile: StabilityProfile,
    regions: List[PromoterRegion],
    outdir: str,
) -> str:
    """一条序列一张独立 PNG"""
    seq_id = profile.seq_id
    short_id = seq_id.split("PM0")[-1][:50] if "PM0" in seq_id else seq_id[:50]

    # 过滤越界预测 + Dmax-ΔG 一致性校验
    valid_regions = []
    mismatch_regions = []
    for r in regions:
        if not r.overlaps_profile(profile):
            logging.warning(f"{seq_id}: promoter {r.start}..{r.end} outside sequence range, skipped")
            continue
        if not r.check_dg_agreement(profile):
            mismatch_regions.append(r)
        valid_regions.append(r)

    fig, ax = plt.subplots(figsize=(16, 5))

    # 自由能折线
    ax.plot(profile.positions, profile.energies,
            color=COLORS["stability_baseline"], linewidth=1.2, alpha=0.9)

    # 均值虚线
    ax.axhline(y=profile.mean_dg,
               color=COLORS["stability_unstable"], linestyle="--",
               linewidth=0.9, alpha=0.7)

    # 均值以下淡色填充（不稳定区的视觉暗示）
    ax.fill_between(profile.positions, min(profile.energies), profile.mean_dg,
                    color=COLORS["stability_unstable"], alpha=0.05)

    # 高亮预测启动子
    for r in valid_regions:
        ax.axvspan(r.start, r.end,
                   alpha=0.22, facecolor=promoter_fill_color(r.dmax),
                   edgecolor=COLORS["border_promoter"],
                   linewidth=0.7, linestyle="--")

        mid = (r.start + r.end) / 2
        y_pos = profile.mean_dg + profile.dg_range * 0.08
        ax.annotate(f'{r.confidence_label()} (Dmax={r.dmax:.1f})',
                   xy=(mid, y_pos), ha="center",
                   fontsize=7, fontweight="bold",
                   color=r.confidence_color(),
                   bbox=dict(boxstyle="round,pad=0.3",
                            facecolor="white", alpha=0.92,
                            edgecolor="#CCCCCC"))

    # Dmax-ΔG 不一致警告
    for r in mismatch_regions:
        mid = (r.start + r.end) / 2
        y_pos = max(profile.energies) + profile.dg_range * 0.02
        ax.annotate("⚠ Dmax/ΔG mismatch",
                   xy=(mid, y_pos), ha="center",
                   fontsize=6, color=COLORS["mismatch_warning"])

    # 统计提示
    if valid_regions:
        hi = sum(1 for r in valid_regions if r.confidence_label() == "High")
        med = sum(1 for r in valid_regions if r.confidence_label() == "Medium")
        low = sum(1 for r in valid_regions if r.confidence_label() == "Low")
        ax.annotate(f'{len(valid_regions)} promoter(s): {hi}H / {med}M / {low}L',
                   xy=(0.98, 0.95), xycoords="axes fraction",
                   ha="right", va="top", fontsize=8, fontstyle="italic",
                   color="#666666")

    # 轴语义标注
    ax.text(0.01, 0.02, "← more stable | less stable →",
            transform=ax.transAxes, fontsize=7,
            verticalalignment="bottom", color="gray")

    ax.set_title(short_id, fontsize=11, fontweight="bold", pad=8)
    ax.set_xlabel("Position (nt)", fontsize=10)
    ax.set_ylabel("Average ΔG (kcal/mol)", fontsize=10)
    ax.legend(["ΔG profile", f"Mean ΔG = {profile.mean_dg:.1f}"],
             loc="upper right", fontsize=9, framealpha=0.85)
    ax.grid(True, alpha=0.25, linestyle=":")

    plt.tight_layout()

    safe_name = seq_id.replace("/", "_").replace("\\", "_").replace(":", "_")[:60]
    outfile = os.path.join(outdir, f"{safe_name}.png")
    plt.savefig(outfile, dpi=180, bbox_inches="tight", facecolor="white")
    plt.close()
    return outfile


def plot_confidence_summary(
    promoters: Dict[str, List[PromoterRegion]],
    outdir: str,
    prefix: str,
) -> str:
    """所有预测按 Dmax 降序的横向柱状图"""
    all_data = []
    for seq_id, regions in promoters.items():
        short_id = seq_id.split("PM0")[-1][:30] if "PM0" in seq_id else seq_id[:30]
        for r in regions:
            all_data.append((short_id, f"{r.start}-{r.end}", r))

    if not all_data:
        logging.warning("No promoter data to summarize.")
        return ""

    all_data.sort(key=lambda x: x[2].dmax, reverse=True)

    labels = [f"{d[0]}\n{d[1]}" for d in all_data]
    dmax_vals = [d[2].dmax for d in all_data]
    colors = [r.confidence_color() for _, _, r in all_data]

    fig, ax = plt.subplots(figsize=(12, max(5, len(all_data) * 0.4)))
    ax.barh(range(len(all_data)), dmax_vals, color=colors,
            height=0.7, edgecolor="white")

    ax.axvline(x=THRESHOLDS["dmax_high"], color=COLORS["confidence_high"],
              linestyle="--", linewidth=1.2, label=f'High (≥{THRESHOLDS["dmax_high"]})')
    ax.axvline(x=THRESHOLDS["dmax_medium"], color=COLORS["confidence_low"],
              linestyle=":", linewidth=1.2, label=f'Low (<{THRESHOLDS["dmax_medium"]})')

    ax.set_yticks(range(len(all_data)))
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_xlabel("Dmax (stability difference)", fontsize=11)
    ax.set_title(f"Promoter Prediction Confidence ({len(all_data)} total)",
                fontsize=13, fontweight="bold")
    ax.legend(fontsize=9)
    ax.invert_yaxis()
    ax.grid(True, alpha=0.3, axis="x")

    plt.tight_layout()
    outfile = os.path.join(outdir, f"{prefix}_confidence_summary.png")
    plt.savefig(outfile, dpi=200, bbox_inches="tight", facecolor="white")
    plt.close()
    logging.info(f"Summary saved: {outfile}")
    return outfile


# ====== 5. 主入口 ======

def main():
    parser = argparse.ArgumentParser(
        description="Visualize PromPredict free energy profiles and promoter predictions",
        epilog="Example: python3 plot_prompredict.py sample_seq",
    )
    parser.add_argument("prefix", help="Prefix of _stb.txt and _PPde.txt files")
    args = parser.parse_args()

    prefix = args.prefix
    stb_file = f"{prefix}_stb.txt"
    ppde_file = f"{prefix}_PPde.txt"

    # 检查文件存在性
    for f in [stb_file, ppde_file]:
        if not os.path.exists(f):
            logging.error(f"File not found: {f}")
            sys.exit(1)

    # 加载
    profiles = load_stb(stb_file)
    promoters = load_ppde(ppde_file)

    total_regions = sum(len(v) for v in promoters.values())
    logging.info(f"Loaded {len(profiles)} sequences, {total_regions} predicted promoters")

    # 报告：有预测但 profile 缺失的序列
    promo_ids = set(promoters.keys())
    profile_ids = set(profiles.keys())
    orphan_promos = promo_ids - profile_ids
    if orphan_promos:
        logging.warning(
            f"{len(orphan_promos)} sequence(s) have promoter predictions "
            f"but no stability profile. These will not be plotted: "
            f"{', '.join(sorted(list(orphan_promos))[:5])}"
            f"{'...' if len(orphan_promos) > 5 else ''}"
        )

    # 报告：有 profile 但无预测的序列
    no_prediction = profile_ids - promo_ids
    if no_prediction:
        logging.info(
            f"{len(no_prediction)} sequence(s) had no predicted promoters "
            f"(will still generate profile plots)"
        )

    # 输出文件夹
    outdir = f"{prefix}_profiles"
    os.makedirs(outdir, exist_ok=True)

    # 逐序列绘图
    for i, (seq_id, profile) in enumerate(profiles.items()):
        regions = promoters.get(seq_id, [])
        out = plot_single_profile(profile, regions, outdir)
        logging.info(f"[{i+1}/{len(profiles)}] {out}")

    # 汇总图
    plot_confidence_summary(promoters, outdir, prefix)

    logging.info(f"Done! All images saved to: {outdir}/")


if __name__ == "__main__":
    main()
