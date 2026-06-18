# 运行流程指南

## 当前状态

### ✓ 已完成
- 所有脚本语法检查通过
- 所有依赖检查通过
- 实现逻辑验证通过
- 测试框架已建立

### ⚠ 需要修复
- 路径设置问题（部分脚本需要修复路径）
- 缺失数据文件（需要运行Step-01生成）

## 运行流程步骤

### 步骤1: 准备数据文件

确保以下文件存在于 `data/` 目录：
- `sample_cluster_4.csv`
- `Genes_17.txt`
- `Genes_249.txt`

### 步骤2: 运行Step-01 (数据预处理)

**R版本**:
```bash
cd /Users/mx/Study/repositories/ImmuCycReg-framework
Rscript code/Step-01_Processdata.R
```

**Python版本**:
```bash
cd /Users/mx/Study/repositories/ImmuCycReg-framework
python code/Step-01_Processdata.py
```

**输出文件**:
- `results/TCGA-LUAD.csv` (或 .Rdata)
- `results/GTEx-LUAD.csv` (或 .Rdata)
- `results/CNV-LUAD.csv` (或 .Rdata)

### 步骤3: 转换RData为CSV (如果需要)

如果R脚本生成了.Rdata文件，需要转换为CSV：

```bash
Rscript code/convert_rdata_to_csv.R
```

### 步骤4: 运行后续步骤

**Step-02a (NMF rank survey / consensus map)**:
```bash
Rscript code/Step-02a_NMF-rank-survey.R
```

This script generates the NMF rank survey and rank-4 NMF consensus map used to evaluate the subtype number.

**Step-02b (final subtype assignment)**:
```bash
Rscript code/Step-02_NMF-ConsensuClustering.R
```

This script reproduces the final K=4 TCGA-LUAD subtype assignment in Supplementary Data 3 S1 and writes `results/NMF/sample_cluster_4.csv`.

**Step-05 (CNV相关性)**:
```bash
Rscript code/Step-05_CNV-Correlation.R
python code/Step-05_CNV-Correlation.py
```

### 步骤5: 验证结果

```bash
# 验证结果一致性
python code/verify_results.py

# 或使用综合测试
python code/comprehensive_test.py
```

## 路径问题修复

### 问题描述
部分脚本在路径设置上有问题，需要确保：
1. 脚本从项目根目录运行
2. 路径正确指向 `data/` 和 `results/` 目录

### 解决方案

**方法1: 从项目根目录运行**
```bash
cd /Users/mx/Study/repositories/ImmuCycReg-framework
python code/Step-05_CNV-Correlation.py
```

**方法2: 修复脚本路径**
所有脚本应使用以下路径设置：
```python
# 获取项目根目录
script_dir = os.path.dirname(os.path.abspath(__file__))
base_dir = os.path.dirname(script_dir)  # 项目根目录
path_read = os.path.join(base_dir, "data")
path_save = os.path.join(base_dir, "results")
```

## 当前可运行的脚本

### 已验证可运行
- Step-05 (已修复路径问题)
- 所有测试脚本

### 需要修复路径
- Step-01 (需要修复路径)
- Step-02 (需要修复路径)
- 其他步骤 (需要修复路径)

## 快速测试

### 测试脚本语法
```bash
python code/comprehensive_test.py
```

### 测试单个脚本
```bash
# 测试Step-05 (已修复路径)
cd /Users/mx/Study/repositories/ImmuCycReg-framework
python code/Step-05_CNV-Correlation.py
```

## 数据文件检查

### 检查现有数据文件
```bash
ls -la data/*.csv data/*.txt
ls -la results/*.csv results/*.Rdata
```

### 检查缺失文件
```bash
# 检查Step-05需要的文件
ls -la results/TCGA-LUAD.csv results/CNV-LUAD.csv
```

## 下一步

1. **修复路径问题**: 确保所有脚本使用正确的路径设置
2. **运行Step-01**: 生成必要的数据文件
3. **运行完整流程**: 按顺序运行所有步骤
4. **验证结果**: 对比R和Python输出

## 注意事项

1. **运行顺序**: 必须按步骤顺序运行（Step-01 → Step-02 → ...）
2. **数据依赖**: 每个步骤依赖前一步的输出
3. **路径要求**: 必须从项目根目录运行脚本
4. **文件格式**: R输出为.Rdata，Python输出为.csv，需要转换对比

## 故障排除

### 问题1: 文件未找到
**错误**: `FileNotFoundError: No such file or directory`
**解决**: 
- 检查文件是否存在
- 确认路径设置正确
- 从项目根目录运行脚本

### 问题2: 路径错误
**错误**: `找不到 data/ 或 results/ 目录`
**解决**: 
- 确保从项目根目录运行
- 检查脚本中的路径设置
- 使用绝对路径

### 问题3: 依赖缺失
**错误**: `ModuleNotFoundError` 或 `Package not found`
**解决**: 
- 安装缺失的Python包: `pip install <package>`
- 安装缺失的R包: `install.packages("<package>")`
