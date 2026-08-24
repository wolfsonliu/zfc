# FASTQ 计数脚本 Python 重写（zfc-count）设计

日期：2026-08-24
状态：已批准（待实现）

## 问题陈述

现有 `bin/library_count_sgrna` 与 `bin/library_count_sgrna_with_barcode` 是
bash + awk 管道实现：每个 read pair 要经过 `cat`×2、`rev|tr`（整文件反向
互补）、`paste`、多段 `awk` 与 `sort`，I/O 和进程开销大；`cat` 读不了
`.gz`；逻辑散落在管道里难以测试和维护。

目标：用 Python 单遍流式重写，原生支持 gzip，单文件跑得更快，输出仍满足
`zfc` 输入契约（`gene\tguide\tbarcode\tcount`，见 `docs/data-format.md`）。

非目标：定长快路径（本文件「后续优化」中的方案 B）、把计数直接接入 `zfc`
主分析流程——仍然保持「计数脚本产出 count 表 → 用户 pivot 成 ctrl/exp →
喂给 `zfc`」的现状。

## 总体方案

新增一个 Python 模块 + 一个统一 CLI，把两个 bash 脚本合并为一种命令的两种
模式（guide / barcode）。单遍、配对流式读取两个 FASTQ，按记录边界分块，
用 `ProcessPoolExecutor` 并行做正则提取 + 计数，主进程合并；library 字典
广播给 worker 实现提前过滤与映射。gzip 用标准库 `gzip`（按 magic bytes
`1f 8b` 自动识别，不依赖扩展名），`python-isal` 加速列为后续可选优化。

## 详细设计

### 文件与接口

- 新模块 `src/zfc/fastqcount.py`：解析、匹配、计数、库映射、`main()`。
- 在 `pyproject.toml` 的 `[project.scripts]` 注册
  `zfc-count = "zfc.fastqcount:main"`（`uv run zfc-count`）。
- v1 仅用标准库（`gzip`/`re`/`argparse`/`collections`/`multiprocessing`），
  不改动 `[project].dependencies`；仅新增 `[project.scripts]` 条目。
- 两个旧 bash 脚本暂时保留，在 README/docs 标注「已被 `zfc-count` 取代
  (deprecated)」；验证稳定后由维护者决定删除。

### CLI

```
zfc-count --mode {guide,barcode} -1/--read1 FILE -2/--read2 FILE \
  [-l/--library FILE] \
  [-p/--guide-pattern REGEX] [-a/--anchor-pattern REGEX] \
  [-b/--barcode-pattern REGEX] \
  [-o/--output FILE|-] [--workers N] [--chunk-size N] [--quiet]
```

- `--mode guide`：只用 `--guide-pattern`（默认 `ACCG([ATGC]{20})GTTT`），
  zfc-ready 模式下 barcode 列输出 `= guide`（对齐 zfc「无 iBAR 时
  barcode=guide」约定）。
- `--mode barcode`：只用 `--anchor-pattern` / `--barcode-pattern`
  （`--guide-pattern` 在本模式不使用，与原脚本一致），默认值与两个 bash
  脚本完全一致：
  - `-a` anchor（sgRNA+barcode 全长）
    `ACCG([ATGC]{20})GTTT[ATGC]{5,35}TGGA([ATCG]{6})AACA`（捕获组 1 = sgRNA
    20nt，捕获组 2 = barcode 6nt）
  - `-p` sgRNA `ACCG([ATGC]{20})GTTT`
  - `-b` barcode `TGGA([ATCG]{6})AACA`
- `-1/-2/-l` 支持 `.gz`（按魔数识别）；输出恒为文本，默认 stdout。
- `-1/-2` 为必填（缺则报错退出），`-l` 为可选（见下）。
- `--workers` 默认 = `min(os.cpu_count() or 1, 8)`；`--workers 1` 用于串行/
  确定性调试。

### 模式：zfc-ready 与 raw

`-l` 是否提供决定输出：

- **zfc-ready 模式（提供 `-l`）**：输出
  `gene\tguide\tbarcode\tcount`，库内条目全量输出、未见到的 key 计 0；
  可直接 pivot 成 ctrl/exp 喂给 `zfc`。
- **raw 模式（省略 `-l`）**：只做提取 + 计数，不做基因映射、不过滤库，
  供开发/QC/提取验证。输出：
  - guide 模式头行 `guide\tcount`
  - barcode 模式头行 `guide\tbarcode\tcount`
  - 文档标注「raw 模式不是 zfc 直接输入（缺 gene 列）」。

library 文件格式（沿用原脚本，见 `docs/data-format.md`）：

- guide 模式：`gene <tab> guide`
- barcode 模式：`gene <tab> guide <tab> barcode`

### 匹配与计数语义（修成更合理语义，逐条文档化）

对每个 read pair `(R1, R2)`，候选位置取 `R1, R2, RC(R1), RC(R2)` 四个
（`RC` = 反向互补；等价于原脚本 `rev|tr` 的双向做法，但不整文件落盘）。

- **guide 模式（zfc-ready）**：按固定优先级
  `R1 → R2 → RC(R1) → RC(R2)` 取第一个命中，捕获组 1 作为 guide，该 pair
  计 1 次（不按双端重复计）。仅统计 library 里存在的 guide；未命中跳过。
- **guide 模式（raw）**：同上提取，但累积所有提取到的 guide（不过滤）。
- **barcode 模式**：在 `R1/R2/RC(R1)/RC(R2)` 用 `--anchor-pattern` 找
  「guide + 同链 barcode」得到 `(G, B1)`，用 `--barcode-pattern` 在对端 read
  找 barcode 得到 `B2`；取库中存在的 `(G, B)`，优先 `B1` 再 `B2`；都不在库
  （zfc-ready）则跳过。每个 pair 计 1 次。raw 模式累积所有提取到的
  `(guide, barcode)`。
- **计数与映射（zfc-ready）**：worker 内用广播的 library 做反向映射
  （`guide→gene` 或 `(guide,barcode)→gene`），只累积库内键；主进程合并各部
  分计数，最后输出所有 library 行（未见到的 key 计 0）。
- **明确修正项**：不再依赖 `library_count_sgrna` 里 `$0 ~ HPATTERN` 的未定义
  变量（原为恒真、靠下游过滤）；方向优先级显式化，替代原脚本隐式的「$1 优先」。

### 性能设计

- 单遍流式：一次读两个文件、按 `4 行/记录` 切块（`--chunk-size`，默认约
  5 万 read-pair），不改写整文件、不重复 pass。
- 并行：生产者按记录边界读块 → `ProcessPoolExecutor.map` → worker 返回部分
  `dict[key]->count`（key 数量 ≈ library 规模）→ 主进程求和。
- gzip：标准库透明解压；`python-isal` 列为后续可选加速。
- 预期：7.5GB / 约 5000 万 read-pair 在 8 核上，主要受 gzip 解压（单线程）+
  正则吞吐约束；会写 benchmark 脚本量化 `reads/s` 与时延。

### 输出

- zfc-ready：头行 `gene\tguide\tbarcode\tcount`；行按 `(gene, guide, barcode)`
  稳定排序（替代 `sort -V`，文档注明排序规则）；仅库内条目，未见 count=0。
- raw：按上述 raw 头行与排序。

### 测试与基准

- `tests/test_fastqcount.py`：合成确定性 FASTQ + library 夹具，覆盖：
  guide/barcode 两模式、raw/zfc-ready 两模式、`.gz` 输入、guide 在 R1/R2 两
  方向、RC 方向、未命中、库外 guide 跳过、库内 count=0、header/排序、配对
  长度不一致报错。
- 基准脚本（非 pytest）跑
  `tmp/sra_tools/fastq/SRR7975595_pass_{1,2}.fastq.gz` 报吞吐。
- 与旧 bash 只在「正常输入下计数数量级可比」层面抽查，不做逐字节 parity
  （因语义修正）。

## 成功标准

- `uvx ruff check`、`uv run pytest` 全绿（新增测试覆盖上述用例）。
- 在 SRR7975595 数据上 `zfc-count` 端到端跑通并输出 headered 表；benchmark
  证明吞吐 ≥ 旧 bash 管线（记录基线数字）。
- 输出能被 `zfc` 现有 `read_raw_count` 消费（无格式破坏）。

## 范围内 / 范围外

- 在范围：guide/barcode 两种计数模式、gzip、并行、统一 CLI、raw/zfc-ready
  两模式、测试与基准、文档更新。
- 不在范围：定长快路径（后续优化·方案 B）、`python-isal`/`pyfastx` 集成、
  把计数接入 `zfc` 主分析流程。

## 开放问题

1. SRR7975595 的 library 尚未提供；开发时用 raw 模式验证提取与吞吐，端到端
   zfc-ready 验证需 library（或合成测试 library）。
2. 是否引入 `python-isal`：v1 先用标准库 gzip（零新依赖），benchmark 后若
   解压占主导再加。
3. 旧 bash 脚本去留：暂时保留 + 标注 deprecated，验证后由维护者决定删除。

## 后续优化（不在本期）

- 定长捕获组（20/6nt）+ 固定侧翼的「直切片哈希计数」快路径，避开正则。
- `python-isal` / `pyfastx` 解压/解析加速。