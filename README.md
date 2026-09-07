# Explicit dissociated sets from tensor lattices

**Erdős 子集和互异问题（Erdős Problem #1）的反证的有效化：显式构造、原始化定理与计算认证。**

> A set of positive integers is *dissociated* if all its subset sums are distinct. Erdős's 1931 conjecture
> `max A >> 2^|A|` was disproved in September 2026 (GPT‑6 Astra; exposition by T. Bloom on erdosproblems.com).
> The disproof is non‑effective because of a "primitive lattice approximation" step. This repository makes that
> step explicit, proves structural theorems about the resulting lattices, and produces computer‑certified explicit
> dissociated sets that beat the previous explicit record (Bohman 1998, `max A ≤ 0.22002·2^n`).

论文（E-JC 投稿格式，单文件）：[`paper/ejc/ejc-main.pdf`](paper/ejc/ejc-main.pdf)（20 页；源码 `paper/ejc/ejc-main.tex`，样式 `e-jc.sty`；用 `tools/bin/tectonic.exe -X compile ejc-main.tex` 编译，无 overfull hbox、无未定义引用、纯 ASCII 源码）。2026‑09‑06 按 PaperSpine 流程重写（引言记分板补全、单调性引理、全部 n 的推论、Table 1 解读、14 条经 Crossref 核对的参考文献），流程产物在 `paper_rewriting_output/`；重写前的稿子保留为 `paper/ejc/ejc-main-draft-2026-09-05.tex/.pdf`。

**目标期刊：Electronic Journal of Combinatorics（E‑JC）**。理由与投稿前必须完成的事项见第 7 节。

---

## 1. 问题与背景

记 `f(n) = inf { max A / 2^(n−1) : A 互异, |A| = n }`。

| 结果 | 界 |
|---|---|
| 2 的幂 | f ≤ 1 |
| Conway–Guy (1968) | N < 0.23513·2^n（n ≥ 40；Bohman 1996 证明整个序列互异） |
| Lunnon (1988) | N < 0.22096·2^n（n ≥ 67，计算机搜索） |
| Bohman (1998) | f ≤ 0.44004（即 N ≤ 0.22002·2^n），本文之前最好的显式构造 |
| Erdős–Moser / Elkies–Gleason / Dubroff–Fox–Xu | f ≫ n^{−1/2} |
| **GPT‑6 Astra (2026‑09)** | 对任意 ε>0，有无穷多 n 使 f(n) < ε（**非有效**） |

Bloom 的阐述指出：非有效性只来自"用原始格逼近 QΛ"这一步，并猜测有效版本应给出 `f(n) ≤ n^{−c/log log n}`。

## 2. 本项目的成果（如实分级）

### 2.1 已完整证明的定理（无条件）

（编号为 E‑JC 版 `ejc-main.pdf` 中的全局编号。）

1. **格的闭式**（Prop. 8）：在 `Z[x_1..x_s]/(x_i^b − 1)` 中，
   `Λ_s = ⊕_{k=1}^{s} 2^{−k}(2+x_1)⋯(2+x_k)(1−x_k)·Z[x_k,…,x_s]`，`v_s = 2^{−s}∏(2+x_k)`。

2. **显式原始化与法向量的张量结构**（Thm. 16）：把单项式 `x_1^{b−1}⋯x_{k−1}^{b−1}x_k^i m` 加到 `Q·Λ_s` 的对应基向量上得到格 `L_s(Q)`，其法向量是
   `u = ũ_1 ⊗ ⋯ ⊗ ũ_s`，`ũ_k ∝ x^{b−1}·adj(φ̄_k)`，层多项式 `φ_k = θ_k(2+x)(1−x) + ρ_k` 由显式递归给出。

3. **精确指标公式**（Thm. 17）：`[sat L : L] = ∏_k cont(adj φ̄_k)^{b^{s−k}}`。于是 d = b^s 维格的原始性证书只需 s 次 b×b 伴随矩阵计算。

4. **层判据与坏素数分类**（Prop. 15, 18）：一层原始 ⟺ 对所有素数 p，`deg gcd_p(φ_k, x^b − 1) ≤ 1`；坏素数只能整除一个仅依赖 b 的固定整数 `b(2^b+1)n_b`，或整除 `gcd(θ_k, ρ_k)`。

5. **第一层的有效定理**（Thm. 19，引言中的 Thm. 3）：3 ∤ b，`rad(P_b) | q` ⟹ `L_1(2q)` 原始，给出显式的 (2q−K_1)-互异整数向量与显式互异集合；Example 20 给出 b=5, Q=16 的完整数值。

6. **显式集合定理**（Thm. 21）：证书 `cont = 1`（全部层）成立时，`u` 的坐标构成 (Q−K_s)-互异集合，且
   `f(ld) ≤ max u / 2^{l(d−1)} ≤ Δ_{b,s}(1 + 10K_s/Q)(Q/2^l)^{d−1}`，其中 `Δ_{b,s} = (1+2^{−b})^{(b^s−1)/(b−1)}(2/3)^s`。

7. **单调性与全部 n**（Lemma 22 + Cor. 23）：经典倍增技巧 `2A ∪ {1}` 给出 `f(n+1) ≤ f(n)`（Conway–Guy、Bohman 都用过；本文只是把它写进 f 的归一化）。因此 Table 1 的每一行都是对所有更大 n 的显式上界：**对每个 n ≥ 43 923 有显式互异集合 max A < 0.158523·2^n；对每个 n ≥ 3 841 966 有 max A < 0.103036·2^n，即 f(n) ≤ 0.206072**。顺带地，反证 + 单调性 ⟹ f(n) → 0 沿所有 n（非有效）。

8. **无穷多好参数**（Thm. 25 + Cor. 26, 27）：在有限可判定条件 H(b,s)（各层多项式 a_k(t), c_k(t) 互素，Def. 24）下，好的 Q 含一个无穷等差数列，因此 `limsup_l f(l·b^s) ≤ Δ_{b,s}`。H(b,s) 已对
   (5,2),(7,2),(11,2),(13,2),(5,3),(7,3),(11,3),(13,3) 认证 ⟹ 无条件地 `lim f(n) ≤ Δ_{13,3} < 0.30299`（仅用 H 证书、不用大伴随计算；显式集合给出更强的 0.206072）。

### 2.2 计算认证的显式构造（精确整数算术）

| b | s | d=b^s | n = l·d | f(n) ≤ | max A / 2^n ≤ | Δ_{b,s} |
|---|---|---|---|---|---|---|
| 7 | 2 | 49 | 1 078 | 0.474515 | 0.237258 | 0.47299 |
| 11 | 2 | 121 | 3 025 | 0.448812 | 0.224406 | 0.44706 |
| **11** | **3** | **1 331** | **43 923** | **0.317045** | **0.158523** | 0.31617 |
| 13 | 3 | 2 197 | 74 698 | 0.304338 | 0.152169 | 0.30299 |
| 17 | 3 | 4 913 | 181 781 | 0.298073 | 0.149037 | 0.29699 |
| 13 | 4 | 28 561 | 1 199 562 | 0.265069 | 0.132535 | 0.26412 |
| **17** | **4** | **83 521** | **3 841 966** | **0.206072** | **0.103036** | 0.20556 |

（完整表见 `results/summary.md` 与论文 Table 1。）这些是据我们所知**第一批显式给出、且超过 Bohman 常数 0.22002 的互异集合**。
n ≤ 25 的小例子全部通过暴力（meet‑in‑the‑middle）验证互异。

### 2.3 条件性结论与猜想

- **猜想 H**（Conj. 28）：H(b,s) 对所有 3 ∤ b 的奇数 b 与所有 s 成立。在此猜想下（Thm. 29，证明只给了概要），对无穷多 n 有 `f(n) ≤ n^{−c/log log n}`，c = 1/20 可取；由于相邻尺寸 n_s 之间只差次多项式因子，单调性把它传到**所有充分大的 n**（Cor. 30，c′ = 1/50）；若好参数具有正密度，则 c = log(3/2) − o(1)。

### 2.4 未解决 / 明确指出的问题

- **构造本身只覆盖稀疏尺寸。** 我们的 n 形如 l·b^s。论文 §7.1 说明三种在构造内部"填充"的手段（追加元素、混合倍增长度、非 2 幂提升）都会损失关于补充元素个数指数级的因子；空隙由经典的倍增引理（Lemma 22）无损填补，"f(n) → 0 是否沿所有 n"不是开放问题（旧稿这一点是错的，已改正）。
- **改进指数。** 指数常数 c = log(3/2) 等于基格高度 h = 3/2 的对数。§7.2（Prop. 32）证明：维数 b 上高度 h ≥ b^α、covol ≤ 1+b^{−Cs} 的可容许对会给出 `f(n) ≤ n^{−α+o(1)}`，且 Vaaler 定理限制 α ≤ 1/2。我们提出确定 `h*(b)` 的问题（Problem 33）。
- **Thm. 29 的完整证明。** 速率定理的伴随/结式高度估计只给了概要；c 的数值依赖于它们，c > 0 的存在性不依赖。

## 3. 诚实的定位

- 反证本身属于 GPT‑6 Astra（经 Bloom 阐述）；本项目提供的是**有效化**：闭式、原始化定理、精确指标公式、有限坏素数集、可认证的显式构造与新的显式记录。
- 这是一篇扎实的后续研究笔记（arXiv 级别，可投中等偏上的数论/组合期刊），**不是**四大期刊级别的成果。四大级别的结果需要例如：无条件的显式衰减速率、多项式衰减、或匹配 n^{−1/2} 的下界。（f(n) → 0 沿所有 n 已由反证 + 倍增引理得到，但那是非有效的。）
- 所有数值断言都可由 `experiments/` 复现；核心证明在论文中完整给出，条件性结论明确标注。

## 4. 复现

环境：Python ≥ 3.12，`numpy`，`sympy`（`pip install -r requirements.txt`）。无需其他依赖；LaTeX 用 `tools/bin/tectonic.exe`（已下载，可选）。

```
python experiments/exp01_base_checks.py      # 格不变量、基对原始化、小例子暴力验证
python experiments/exp02_gcd_patterns.py     # gcd 模式（3|b 的系统性失败等）
python experiments/exp03_tensor_structure.py # 可容许性抽检、多项式判据 vs 直接 gcd
python experiments/exp04_tensor_check.py     # 张量递归 vs 直接子式（方向、指标公式）
python experiments/exp05_polynomials.py      # 层多项式、a_1 = N − (q/b)N'、结式
python experiments/exp06_certificates.py     # 定理A验证、H(b,s) 认证、好Q搜索
python experiments/exp07_K_exact.py          # K_s 递归 vs 精确对偶基范数
python experiments/exp08_final_checks.py     # s=2 暴力验证、真实阈值、(17,3),(13,4),(17,4)
python experiments/make_certificates.py      # 生成 results/certificates.json, results/summary.md（加 --big 含 (17,4)）
```

核心代码：`src/dissociated.py`（互异性检验）、`src/lattice.py`（可容许对、子式、对偶基）、`src/tensor.py`（多项式模型、层递归、纯整数 Bareiss 伴随、K_s、Δ）。

## 5. 目录

```
paper/ejc/    ejc-main.tex（单文件，E-JC 格式）, ejc-main.pdf, e-jc.sty, ejc-sample.tex, ejc-main-draft-2026-09-05.*（重写前的稿子）
paper_rewriting_output/  PaperSpine 流程产物（研究档案、引文库、写作理由矩阵、审稿人审计、final_paper/main.tex 与 PDF）
src/          核心库
experiments/  可复现实验（exp01–exp08, make_certificates）
results/      认证表与实验输出
literature/   erdosproblems.com 与 E-JC 投稿指南的纯文本存档
tools/        抓取/汇总脚本，tectonic 可执行文件
```

## 7. 投稿：目标期刊与待办

**首选：Electronic Journal of Combinatorics (E‑JC)。**
- 主题延续性：Bohman 1998 年的构造（本文超越的显式记录）就发表在 E‑JC 5 (1998) R3；互异集合/Erdős 问题是 E‑JC 的经典题材。
- 明确接受计算机辅助证明，无篇幅限制，开放获取，中位审稿周期约 6 个月，无版面费。
- 有成文的 **AI 政策**（`literature/raw/ejc_about.txt`）：允许 AI 协助，但作者对全部内容负责，必须亲自核对所有证明与计算，保证引用真实。

**备选**：Journal of Number Theory（几何数论/有效结果的常规去处）；Integers（组合数论，审稿快）；Research in Number Theory；Experimental Mathematics（若想突出计算与猜想 H 的实验证据）。

**投稿前必须由作者本人完成：**
1. 逐行核对论文中的每个证明（尤其 Thm. 17 指标公式、Prop. 18 坏素数分类、Thm. 25 的归纳、新增的 Lemma 22 / Cor. 23 / Cor. 30），并重跑 `experiments/` 复现所有数字；最好用另一套系统（PARI/GP 或 Magma）独立复算一个证书。
2. 参考文献的卷期页码与 DOI 已于 2026‑09‑06 经 Crossref API 核对（14 条，见 `paper_rewriting_output/citation_support_bank.md`）；投稿前再确认 Bloom 页面是否已有稳定版本（arXiv）。
3. 占位符：作者、单位、邮箱与致谢中的 AI 使用声明已填（2026‑09‑07）；仍需填写代码仓库 URL（建议 GitHub + Zenodo DOI）；`\dateline` 由投稿/录用日期决定。致谢中"the author has checked all proofs and computations"一句必须在投稿前成为事实（见第 1 条）。
4. 先挂 arXiv（E‑JC 允许），并检索是否已有他人发布有效化版本——这是两天前的热点结果，撞车风险真实存在。
5. E‑JC 初投只需 PDF；录用后上传单一 `.tex` 源文件（已满足：单文件、无 overfull hbox）。

## 6. 参考

- T. F. Bloom, *Erdős Problem #1*, https://www.erdosproblems.com/1（含 GPT‑6 Astra 构造的阐述，2026‑09‑03），访问于 2026‑09‑04。
- T. Bohman, *A construction for sets of integers with distinct subset sums*, Electron. J. Combin. 5 (1998), R3, doi:10.37236/1341.
- T. Bohman, *A sum packing problem of Erdős and the Conway–Guy sequence*, Proc. Amer. Math. Soc. 124 (1996), 3627–3636.
- W. F. Lunnon, *Integer sets with distinct subset-sums*, Math. Comp. 50 (1988), 297–320.
- Q. Dubroff, J. Fox, M. W. Xu, *A note on the Erdős distinct subset sums problem*, SIAM J. Discrete Math. 35 (2021), 322–324.
- S. Steinerberger, *Some remarks on the Erdős distinct subset sums problem*, Int. J. Number Theory 19 (2023), 1783–1800.
- T. Horesh, Y. Karasik, *Equidistribution of primitive lattices in R^n*, Q. J. Math. 74 (2023), 1253–1294.
- I. Aliev, *Siegel's lemma and sum-distinct sets*, Discrete Comput. Geom. 39 (2008), 59–66.
- J. D. Vaaler, *A geometric inequality with applications to linear forms*, Pacific J. Math. 83 (1979), 543–553.
- 完整 14 条见论文参考文献。
