# E-JC 投稿操作清单（2026-09-08 版）

投稿入口：https://www.combinatorics.org/ojs/index.php/eljc/about/submissions （需注册 OJS 账号并登录）。

## 一、初投（只需 PDF）

1. 登录后选择 "Make a new submission"，栏目选研究论文（Research Paper）。
2. 投稿清单逐项勾选（原文见 `literature/raw/ejc_submissions.txt`，不入库）：未在他处发表或审稿中；按作者指南准备、PDF 格式、含摘要；用自己的语言写成并给出全部来源的应有 credit；同意录用后用期刊 LaTeX 样式做最终排版；全部作者同意投稿；**已阅读并遵守 AI 政策**；同意编辑用 AI 做最终检查。
3. 上传 `paper/ejc/ejc-main.pdf`（28 页）。不要在此阶段上传源码。
4. 元数据：
   - 标题：Explicit dissociated sets from tensor lattices: primitivity certificates and an effective decay rate
   - 作者：Zhiyu Liu；单位：Tianjin Medical University General Hospital, Tianjin Medical University, Tianjin, China；邮箱：07eb1f@gmail.com（系统要求直接输入带重音的字符，本单位无此问题）
   - HTML 摘要：粘贴 `paper/ejc/submission_abstract.txt` 中的摘要段（已去掉 `\emph`、自定义宏和 `\,` 千位分隔，只保留 `$...$` 数学）
   - MSC：11B75, 11H06, 11H31, 11Y16
5. 提交后记下投稿编号；E-JC 中位审稿周期约 6 个月。

## 二、AI 政策对应关系

- 致谢中已披露：Claude（Anthropic，经 Cursor）与 Codex（OpenAI）参与了论证的开发与复核、计算的实现与检查、文本的起草与修改；作者对证明、计算与参考文献负全责；新结果未在 Lean 中形式化。
- 政策要求"自己核对并给出足够细节供他人核查"：全部证明完整写出；15 行证书有压缩精确见证与独立验证器（`python tools/verify_certificates.py`），第二条独立路径为 `python tools/second_check_table.py`。

## 三、录用后的最终排版（现在不要做）

1. 删除 `ejc-main.tex` 序言中的这段覆盖代码以及 `\pagestyle{plain}`：
   ```
   \makeatletter
   \def\the@dateline{...}
   \renewcommand{\ps@plain}{...}
   \makeatother
   \pagestyle{plain}
   ```
   改用标准 `\dateline{投稿日期}{录用日期}{TBD}`（格式如 `Jan 1, 2024`），期刊填第三项。
2. 保持 `e-jc.sty` 原样；重新用 tectonic 或 pdflatex 编译，确认无 overfull hbox。
3. 上传单一 `.tex` 源文件（含内置 `thebibliography`）与 PDF 作为补充文件。

## 四、可选：Zenodo DOI

1. 用 GitHub 账号登录 https://zenodo.org ，在 GitHub 集成页面打开仓库 `Chth1z/math-erdos1-dissociated-sets` 的开关。
2. 之后在 GitHub 新建一个 release（例如 `v1.0.1`）；Zenodo 会自动归档并生成 DOI。
3. 把 DOI 写进论文致谢的仓库一句，并在录用后的最终版中体现。

## 五、arXiv（可选，非必需）

E-JC 不要求预印本。若想挂 arXiv：用 tmu.edu.cn 邮箱注册可能获得自动认可；否则需要 math.NT 或 math.CO 的担保人。不建议使用 viXra。
