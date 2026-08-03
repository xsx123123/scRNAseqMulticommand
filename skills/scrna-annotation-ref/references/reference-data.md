# 参考数据获取指引（scrna-annotation-ref）

> 红线：以下所有文件**禁止打进技能包**，一律放仓库 `Celldex/` 目录或共享数据卷。路径约定见仓库根 `scRNAseqMulticommand.yaml` 的 `singeler_reference` 节（节名拼写如此，照抄）。

## 一、SingleR 参考 rds（7 个，仓库不含，需自行下载）

主流程直接 `readRDS()` 这些文件（`src/cli/00.Initialize_Project.r` 按 `-A` 参数加载选定的一个；`src/core/10.annotation.r` 的 `RunSingleR_Unified` 则无条件加载全部 7 个）。文件来源是 Bioconductor **celldex** 包（数据经 ExperimentHub 分发）， celldex 数据函数返回 `SummarizedExperiment`，`saveRDS()` 存成下表文件名即可。

映射关系（celldex 函数名 → 仓库要求的文件名；仓库代码内无 celldex 调用，此映射依据 celldex 官方数据集命名与 yaml 声明的文件名对应）：

| 物种 | yaml 声明的文件名 | celldex 下载函数 | 数据集说明 |
|---|---|---|---|
| 人 | `Celldex/HumanPrimaryCellAtla.rds` | `celldex::HumanPrimaryCellAtlasData()` | HPCA，人原代细胞图谱，通用首选（-A 默认） |
| 人 | `Celldex/HumanBlueprintEncode.rds` | `celldex::BlueprintEncodeData()` | Blueprint/ENCODE 血液 bulk 参考（-A 传错拼 `HuamnBlueprintEncode`） |
| 人 | `Celldex/HumanDICEImmuneCell.rds` | `celldex::DatabaseImmuneCellExpressionData()` | DICE 免疫细胞表达库 |
| 人 | `Celldex/HumanMonacoImmune.rds` | `celldex::MonacoImmuneData()` | Monaco 外周血免疫精细亚型 |
| 人 | `Celldex/HumanNovershternHematopoietic.rds` | `celldex::NovershternHematopoieticData()` | Novershtern 造血干祖细胞 |
| 鼠 | `Celldex/MouseRNA.rds` | `celldex::MouseRNAseqData()` | 小鼠 RNA-seq 参考（-A 默认 `MouseRNAref`） |
| 鼠 | `Celldex/MouseImmGen.rds` | `celldex::ImmGenData()` | ImmGen 小鼠免疫细胞 |

下载示例（需要能访问 ExperimentHub 的网络环境，首次下载会建本地缓存）：

```r
BiocManager::install("celldex")   # 版本要求见 environment.md
library(celldex)
saveRDS(HumanPrimaryCellAtlasData(),        "Celldex/HumanPrimaryCellAtla.rds")
saveRDS(BlueprintEncodeData(),              "Celldex/HumanBlueprintEncode.rds")
saveRDS(DatabaseImmuneCellExpressionData(), "Celldex/HumanDICEImmuneCell.rds")
saveRDS(MonacoImmuneData(),                 "Celldex/HumanMonacoImmune.rds")
saveRDS(NovershternHematopoieticData(),     "Celldex/HumanNovershternHematopoietic.rds")
saveRDS(MouseRNAseqData(),                  "Celldex/MouseRNA.rds")
saveRDS(ImmGenData(),                       "Celldex/MouseImmGen.rds")
```

要点：

- **文件名必须与上表完全一致**（yaml 按文件名读，错了 `readRDS` 直接报错）；
- 注意 `HumanBlueprintEncode.rds` 文件名本身是正确拼写，错拼只出现在 `-A` 参数值与 yaml 键名里；
- 参考对象的 `label.main` / `label.fine` 列必须保留，SingleR 注释用 `label.main`；
- Docker 运行时把准备好的 Celldex 目录挂载到容器内对应路径。

## 二、ScType / marker 数据库文件（4 个，仓库自带）

随仓库分发在 `Celldex/` 下，缺失时用 `git checkout` / 重新 clone 恢复，**不要**从技能包复制：

| 文件 | 大小 | 内容 | 用途 |
|---|---|---|---|
| `Cell_marker_Human.txt` | ~38M | CellMarker 2.0 人 marker 表，TSV 20 列 | `-F Cellmarker` 时的人源 marker 参考（`fread` 读取） |
| `Cell_marker_Mouse.txt` | ~10M | CellMarker 2.0 鼠 marker 表，TSV 20 列 | `-F Cellmarker` 时的鼠源 marker 参考 |
| `PanglaoDB_markers_27_Mar_2020.tsv` | ~1.2M | PanglaoDB marker 表（2020-03-27 快照） | `-F PanglaoDB` 时的 marker 参考（`read.table` 读取） |
| `ScTypeDB_full.xlsx` | ~32K | ScType 内置 cell-type×marker 库（正/负 marker 集） | ScType 打分底库，`gene_sets_prepare()` 按 `-O` 器官过滤 |

- `-O/--organ` 的值必须出现在所用库的 tissue/organ 字段中，否则该组织整轮 ScType 跳过（日志 warning，不报错）；
- ScType 打分阈值：cluster 得分 < cluster 细胞数/4 → 判 `Unknown`（见 `src/core/10.annotation.r`）。

## 三、CellID

CellID 走 Bioconductor **CellID** 包内置的 MCA（Mouse Cell Atlas）/ HCA 参考签名，函数内自动经 ExperimentHub 获取，**无需手工准备文件**，只需要包装好且能访问 ExperimentHub 缓存。

## 四、历史遗留边界

`tools/SingleR.r` 是 2024.8 的历史函数库，参考路径硬编码 `/glusterfs/...`，已被 `src/core/10.annotation.r` 取代，主流程不再 source 它。不要修复它、不要以它为准。
