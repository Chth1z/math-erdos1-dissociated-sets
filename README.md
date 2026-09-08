# Explicit dissociated sets from tensor lattices

**Primitivity certificates and an effective decay rate**

本项目研究 Erdős 子集和互异问题中的显式构造。论文给出保留张量结构的格扰动、精确饱和指标公式、局部原始性判据及可复核的数值构造，并证明无条件的有效定量衰减率。

论文：[E-JC 格式 PDF](paper/ejc/ejc-main.pdf) · [单文件 LaTeX 源码](paper/ejc/ejc-main.tex)。目标期刊为 Electronic Journal of Combinatorics。作者为 Zhiyu Liu；单位与联系方式见论文。当前文件是投稿前修订稿。

## 结果与适用范围

正整数集合 A 的所有子集和互不相同时，称 A 为 dissociated。记

```
f(n) = inf { max A / 2^(n-1) : A dissociated, |A| = n }.
```

本文的无条件结构结果包括：

- 在 `Z[x_1,...,x_s]/(x_i^b-1)` 中给出迭代格的闭式。这里 b≥5 为奇数且 3∤b，d=b^s。
- 按固定规则向缩放格的基向量加入单位向量，得到 `L_s(Q)`；其原始整数法向量分解为 `u = u_1 ⊗ ... ⊗ u_s`。
- 各层向量由 b×b 循环矩阵的伴随计算给出。若 `gamma_k` 为第 k 层伴随系数的 gcd，则
  `[sat L_s(Q) : L_s(Q)] = ∏ gamma_k^(b^(s-k))`。
- 原始性归约为有限域上的多项式 gcd 判据；除了同时整除层参数的素数，坏素数来自只依赖 b 的固定有限集。
- 当 Q 满足显式尺度条件且所有 gamma_k=1 时，u 的坐标是 `(Q-K_s)`-dissociated。二进制提升给出 n=ld 个元素，以及精确有理数界 `f(ld) ≤ max u / 2^(l(d-1))`。

这里的“认证”指精确整数恒等式、gcd 与整数不等式检查；本文的新定理尚未在 Lean 中形式化。

## 已认证的具体上界

以下小数均从精确有理数**向上取整至六位小数**，未先转成浮点数。完整15行数据见 [results/summary.md](results/summary.md)。

| b | s | d | n=ld | f(n) ≤ | max A / 2^n ≤ |
|---:|---:|---:|---:|---:|---:|
| 7 | 2 | 49 | 1 078 | 0.474515 | 0.237258 |
| 11 | 2 | 121 | 3 025 | 0.448813 | 0.224407 |
| 11 | 3 | 1 331 | 43 923 | 0.317045 | 0.158523 |
| 13 | 3 | 2 197 | 74 698 | 0.304339 | 0.152170 |
| 17 | 3 | 4 913 | 181 781 | 0.298073 | 0.149037 |
| 13 | 4 | 28 561 | 1 199 562 | 0.265070 | 0.132535 |
| 17 | 4 | 83 521 | 3 841 966 | 0.206073 | 0.103037 |

这些集合改善 Bohman (1998) 的显式常数 `max A ≤ 0.22002·2^n`。经典操作 `A ↦ 2A ∪ {1}` 保持 `max A / 2^|A|`，因此每行上界适用于所有更大的 n。特别地，对每个 n≥3 841 966，都有本文构造的显式集合满足

```
max A < 0.103037 * 2^n,     f(n) <= 0.206073.
```

2026-09-08 的修复纠正了旧稿的最近舍入错误。旧值 `0.206072` 和 `0.103036` 不能由该行证书支持；这两个更小的小数不应继续作为该构造的上界引用。

## 首层非消失定理与无条件速率

设 `g(x)=2-x-x^2`，`psi=adj(conjugate(theta*g+rho))`，定义 b−1 次齐次型

```
A_b(theta,rho) = psi_0,
C_b(theta,rho) = 2*psi_1 + psi_2,
r_b = Res_q(A_b(q,1), C_b(q,1)).
```

**正文现已证明：对所有奇数 b≥5 且 3∤b，r_b 都非零。** 先前保留的全称猜想已由一般证明解决，不再作为假设。

证明的关键归约是：若存在公共复根，则循环二阶递推给出复数 r,w，满足

```
r+w = -1,
r^(b-1)*(r-1) = w^(b-1)*(w-1) = 1.
```

r 是非零代数整数。对它的每个共轭 z，另一个数仍为 −1−z；比较两条乘积等式的模长，迫使 Re z=−1/2，继而 |z|<1。于是 r 的所有共轭的模都小于1，与其范数为非零整数矛盾。正文还完整处理 q=0、奇异循环矩阵及重复特征根。

齐次递推将首层互素性传到全部深度，所以 **H(b,s) 对所有允许的 b 和全部 s≥1 成立**。取

```
M_b = rad(P_b) * abs(r_b),      log M_b <= 4*b^3,
```

即可使所有满足尺度门槛的 `Q=2^(s+1)t`（M_b 整除 t）给出原始格。这里不依赖 s 的是 t 的模数 M_b；Q 的步长仍是 `2^(s+1)M_b`。

结合正文已证明的模数高度、参数取整和尺寸插值估计，现在**无条件**得到：

- 无穷多个 n 满足 `f(n) ≤ n^(-1/(20 log log n))`；
- 所有充分大的 n 满足 `f(n) ≤ n^(-1/(50 log log n))`。

前面的15行有限证书与参数保持有效。首层结式的四个完整数值，以及更大范围的精确多项式gcd计算，作为一般定理的独立检查。关于更高基格高度是否能带来多项式衰减的问题仍然开放；正文明确列出该方向所需的原始化尺度和尺寸间距前提。

## 与原反证的关系

2026年公开的 GPT-6 Astra 反证已有 Lean 验证；原证明使用 Smith 型整基变换、饱和双对角扰动及具体权重。Bloom 的网页阐述另采用原始格等分布的定性形式。Horesh–Karasik 的等分布论文有定量误差，不能把它本身称作非有效定理。

本文比较的是保留张量坐标的扰动：法向量和饱和指标可由 s 个 b 维计算获得，并有明确的扰动常数与可实际复算的记录。首层非消失定理及统一模数估计给出了该构造的无条件有效衰减率。原反证加经典单调性已经给出所有 n 上的 f(n)→0；极限趋零本身不是本项目的新结果。

原始来源：[Bloom 的阐述](https://www.erdosproblems.com/1)、[原 Lean 仓库固定版本](https://github.com/tadamcz/erdos1/tree/0e395153306f34b3829d118b85bdd704136f2843)。论文参考文献另列 Bohman、Horesh–Karasik、Vaaler 等原始文献。

## 复现

环境：Python 3.12 或更新版本，依赖见 [requirements.txt](requirements.txt)。从仓库根目录运行：

```sh
python -m pip install -r requirements.txt
python experiments/make_certificates.py
python tools/verify_certificates.py
python experiments/exp09_certificate_checks.py
python experiments/exp10_coprimality.py
python tools/tex_static_check.py paper/ejc/ejc-main.tex
```

生成器默认覆盖全部15行，包括 (17,4)。manifest 与精确见证位于 `results/`，存储格式见 [CERTIFICATE_FORMAT.md](results/CERTIFICATE_FORMAT.md)。生成器和独立验证器采用不同的伴随/行列式路线；校验和用于确认文件和整数的一致性，数学认证还会重建整数、检查行列式与伴随恒等式、gcd、层递归、尺度条件和小数上界。

更多实验：

| 脚本 | 检查内容 |
|---|---|
| `experiments/exp01_base_checks.py` | 基格不变量、原始性与小集合检查 |
| `experiments/exp02_gcd_patterns.py` | gcd 模式和坏参数 |
| `experiments/exp03_tensor_structure.py` | 张量不变量、有限参数的可容许性抽检和有限域判据 |
| `experiments/exp04_tensor_check.py` | 直接最大子式与张量法向量、饱和指标的精确比较 |
| `experiments/exp05_polynomials.py` | 层多项式、导数恒等式和结式 |
| `experiments/exp06_certificates.py` | 第一层定理、H 条件的取值见证与参数搜索 |
| `experiments/exp07_K_exact.py` | K_s 的对偶基计算 |
| `experiments/exp08_final_checks.py` | 小集合的穷举检查与较大参数 |
| `experiments/exp09_certificate_checks.py` | 四套算法交叉核对、拒绝篡改证书、首层结式的完整多项式验证 |
| `experiments/exp10_coprimality.py` | 首层非消失定理的精确gcd、递推恒等式、重复根与3整除维数的对照检查 |

有限枚举不是针对所有实参数 t 的证明；基对的 (A1)、(A2) 已在正文给出自包含证明。部分小例子低于定理的充分尺度门槛，它们由独立穷举认证。

LaTeX 使用仓库已有的 Tectonic 或兼容的 LaTeX 环境。PowerShell 中：

```powershell
Set-Location paper/ejc
../../tools/bin/tectonic.exe -X compile ejc-main.tex
```

## 文件与修订

- `paper/ejc/ejc-main.tex`、`ejc-main.pdf`：当前主稿。
- `src/`：格、张量递归和互异性检查。
- `experiments/`、`tools/verify_certificates.py`：生成与独立复核。
- `results/`：精确证书、校验和、表格与实验输出。
- `paper/ejc/ejc-main-draft-2026-09-05.*`：历史草稿，不代表当前结论。

2026-09-08 修订先补全了速率证明，随后证明首层结式对所有允许b非零，消去了该速率的条件假设；同时纠正一般提升与 padding 的归一化公式，补充原证明的固定版本引用，并更新 AI 使用披露。内部旧评审和旧流程日志保留为历史记录，不作为当前稿件或认证数据。
