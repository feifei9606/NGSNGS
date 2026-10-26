#!/usr/bin/env python3
"""把 ART 的 Illumina 质量剖面转换成 NGSNGS 的 AccFreq 质量剖面。

背景
----
NGSNGS 的 -q1/-q2 要求「每行列数恒定」的 AccFreq 格式：

    line0: 质量分集合 Q（升序）
    line1: 对应的错误概率 10^(-q/10)
    随后 5 个区块 × cycle 行，区块顺序固定为 A, T, G, C, N，
    每行的值 = 该碱基在该 cycle 上取各个质量分的**权重**（相对频率）。
    权重语义见 RandSampling.cpp 的 ransampl_set：它对 p[] 做 sum 归一化
    （P[i] = p[i]*n/sum），因此行值不是 CDF，而是各档的相对频率。

ART 的剖面是「每 cycle 两行」的**累积计数**格式：

    <base>  <pos>  q1 q2 ... qk        # 该 cycle 的质量分档（升序）
    <base>  <pos>  c1 c2 ... ck        # 累积计数，c_k = 总数
    <base> ∈ {'.', A, T, G, C, N}，'.' 为不分碱基的总剖面

直接把 ART 表喂给 NGSNGS 会失败（实测报
`Problems, number of columns is different nbins: 7 at: 8` 并产出 0 字节）：
① 各 cycle 的质量档位数不同 → 列数不一致；
② ART 是累积计数、NGSNGS 要权重，语义不同。

NGSNGS 自带的 misc/ReadQualConverter.cpp 正是为此而写，但它在本表上有 bug
（CDF_array[col_i-1] 在 col_i=0 处越界，产出 `inf`；且未处理各 cycle 档位数
与全局档位数不一致、以及本表 A 块第 75 位为空的情况）。故此处重写。

做法
----
- 全局质量分集合 = 所有区块 q 的并集；
- 对每个 (碱基, cycle)：把累积计数差分成各档频率，映射到全局集合
  （不在该 cycle 档位中的质量分权重为 0）；
- 空行（本表 A 块 pos=75）回退到 '.' 区块同 cycle；仍缺则回退到该碱基
  上一个有数据的 cycle。

用法::

    python3 art2ngsngs.py -i NextSeq500v2L75R1.txt -o NextSeq500v2L75R1_NGSNGS.txt
"""

from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

BASES = ["A", "T", "G", "C", "N"]  # NGSNGS 要求的区块顺序


def parse_art(path: Path):
    """读 ART 剖面，返回 {base: {pos: (qs, cumulative_counts)}}。"""
    raw = {}
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            parts = line.split("\t")
            while parts and parts[-1] == "":
                parts.pop()
            if len(parts) < 2:
                continue
            base, pos, vals = parts[0], int(parts[1]), parts[2:]
            raw.setdefault(base, []).append((pos, vals))

    blocks = {}
    for base, rows in raw.items():
        per_pos = {}
        for i in range(0, len(rows) - 1, 2):
            (p1, v1), (p2, v2) = rows[i], rows[i + 1]
            if p1 != p2:
                raise ValueError(f"{path}: 行不成对 base={base} pos {p1}!={p2}")
            if not v1 or not v2:
                continue  # 空行，稍后回退
            qs = [int(x) for x in v1]
            cs = [int(x) for x in v2]
            if len(qs) != len(cs):
                raise ValueError(f"{path}: base={base} pos={p1} 档位数不一致 {len(qs)}!={len(cs)}")
            per_pos[p1] = (qs, cs)
        blocks[base] = per_pos
    return blocks


def bin_freqs(cs):
    """累积计数 -> 各档频率（按顺序）。"""
    total = cs[-1]
    if total <= 0:
        return None
    freqs = [cs[0]] + [cs[i] - cs[i - 1] for i in range(1, len(cs))]
    return [f / total for f in freqs]


def convert(path: Path, out_path: Path) -> dict:
    blocks = parse_art(path)
    missing_bases = [b for b in BASES if b not in blocks]
    if missing_bases:
        raise ValueError(f"{path}: 缺少区块 {missing_bases}")

    alphabet = sorted({q for b in BASES for (qs, _c) in blocks[b].values() for q in qs})
    n_cycles = max(p for b in BASES for p in blocks[b]) + 1

    # 预计算各 (base,cycle) 的权重；带回退
    def weights_for(base: str, cycle: int):
        for cand_base in (base, "."):
            hit = blocks.get(cand_base, {}).get(cycle)
            if hit and hit[1][-1] > 0:
                qs, cs = hit
                w = {q: 0.0 for q in alphabet}
                for q, f in zip(qs, bin_freqs(cs)):
                    w[q] += f
                if sum(w.values()) > 0:
                    return w
        return None

    rows = []
    fallback_used = []
    for base in BASES:
        last_good = None
        for cycle in range(n_cycles):
            w = weights_for(base, cycle)
            if w is None:
                if last_good is None:
                    raise ValueError(f"{path}: base={base} cycle={cycle} 无可用数据且无回退")
                w = last_good
                fallback_used.append((base, cycle, "prev-cycle"))
            else:
                last_good = w
            rows.append(w)

    with open(out_path, "w") as out:
        out.write("\t".join(str(q) for q in alphabet) + "\t\n")
        out.write("\t".join(f"{10 ** (-q / 10.0):.6f}" for q in alphabet) + "\t\n")
        for w in rows:
            out.write("\t".join(f"{w[q]:e}" for q in alphabet) + "\t\n")

    return {
        "alphabet": alphabet,
        "n_cycles": n_cycles,
        "rows": len(rows),
        "blocks": {b: len(blocks[b]) for b in BASES},
        "fallback_used": fallback_used,
    }


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description="ART -> NGSNGS 质量剖面转换")
    p.add_argument("-i", "--input", type=Path, required=True, help="ART 剖面 .txt")
    p.add_argument("-o", "--output", type=Path, required=True, help="NGSNGS AccFreq .txt")
    args = p.parse_args(argv)

    info = convert(args.input, args.output)
    print(f"输入 : {args.input}")
    print(f"输出 : {args.output}")
    print(f"质量分集合 : {info['alphabet']}")
    print(f"cycle 数   : {info['n_cycles']}")
    print(f"输出行数   : {2 + info['rows']}")
    print(f"各区块档位数 : {info['blocks']}")
    if info["fallback_used"]:
        print(f"回退行     : {info['fallback_used']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
