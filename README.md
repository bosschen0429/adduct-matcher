# ESI加合物比對程式使用說明

## 📋 程式功能

此程式用於從質譜數據中識別同一化合物的不同加合物(Adduct)形式，支援23種常見的ESI加合物。

### 核心原理
1. **滯留時間(RT)分組**: 找出在相近滯留時間內的所有訊號
2. **質量差比對**: 計算不同m/z值之間的質量差（大減小）
3. **加合物識別**: 將質量差與已知的加合物質量差異表比對
   - PPM容許範圍 = ppm × 較大的m/z值 / 1,000,000
   - 例如：20 ppm 對 m/z 264 = ±0.00528 Da
4. **多重配對處理（方案2）**: 
   - **Base有多個Pair候選** → 選擇RT最接近的
   - **Pair有多個Base候選** → 保留所有配對，列出所有可能性
   - 在Adduct_Type和Matched_Adduct欄位中用分號`;`分隔
   - 按PPM誤差排序，最可能的配對列在前面
5. **結果註記**: 標示出可能的加合物類型並保留所有原始欄位

---

## ✨ 主要特色

- ✅ **自動檔案讀取**: 支援 Excel (.xlsx, .xls, .xlsm, .xlsb), CSV, TSV
- ✅ **智能欄位識別**: 自動識別 RT, m/z, Intensity 欄位（支援多種名稱變體）
- ✅ **保留所有欄位**: 原始數據的所有欄位都會被保留在結果中
- ✅ **完整Excel報告**: 自動生成包含4個工作表的分析報告
- ✅ **23種ESI加合物**: 涵蓋最常見的正離子模式加合物
- ✅ **精確比對**: 可調整PPM和RT容許範圍

---

## 🚀 快速開始

### 安裝需求

```bash
pip install pandas openpyxl numpy --break-system-packages
```

### 方法1: 使用圖形介面（推薦）

**最簡單的使用方式！**

1. 執行程式:
```bash
python adduct_matcher.py
```

2. 在彈出的視窗中:
   - 點擊「選擇檔案」按鈕
   - 選擇你的質譜數據檔案（Excel, CSV, 或 TSV）
   - 設定參數（或使用預設值）
   - 點擊「開始比對」

3. 完成！結果會自動儲存為 `檔名_adduct_results.xlsx`

### 方法2: 在程式中使用

```python
from adduct_matcher import AdductMatcher

# 建立比對器（預設20 ppm容許誤差）
matcher = AdductMatcher(ppm_tolerance=20.0)

# 處理檔案（一行搞定！）
results = matcher.process('你的檔案.xlsx', rt_tolerance=0.05)
```

---

## 📁 支援的檔案格式

### Excel格式
- `.xlsx` (Excel 2007+)
- `.xls` (Excel 97-2003)
- `.xlsm` (Excel 啟用巨集)
- `.xlsb` (Excel 二進位)

### 文字格式
- `.csv` (逗號分隔)
- `.tsv` (Tab分隔)
- `.txt` (Tab分隔)

---

## 🔍 自動欄位識別

程式會自動識別以下欄位（不區分大小寫）:

### RT (滯留時間) - **必需**
可識別的欄位名稱:
```
rt, RT, RT (min), RT(min)
retention time, Retention Time, Retention Time (min)
time, Time
r.t., R.T., R.T.
保留時間, 滯留時間
```

### m/z (質荷比) - **必需**
可識別的欄位名稱:
```
m/z, M/Z, mz, MZ
mass, Mass
mass to charge, Mass to Charge
precursor ion m/z, Precursor m/z
質量, 質荷比
```

### Intensity (強度) - **必需**
可識別的欄位名稱:
```
intensity, Intensity, Int
abundance, Abundance, Abund
height, Height
area, Area
precursor ion intensity, Precursor Intensity
強度, 峰高, 峰面積
```

**其他欄位**: 所有其他欄位都會被自動保留在結果中

---

## 📋 輸入檔案格式

### 最簡單的格式
| RT (min) | m/z      | Intensity |
|----------|----------|-----------|
| 10.99    | 242.1191 | 50000     |
| 10.99    | 264.1010 | 38000     |

### 包含額外欄位的格式（推薦）
| Peak ID | RT (min) | m/z      | Intensity | Area  | S/N | 其他欄位... |
|---------|----------|----------|-----------|-------|-----|------------|
| 10      | 10.99    | 242.1191 | 50000     | 250k  | 200 | ...        |
| 15      | 10.99    | 264.1010 | 38000     | 190k  | 150 | ...        |

**重點**: 
- 必須包含RT、m/z、Intensity三個欄位
- 欄位名稱會自動識別，不用擔心大小寫或格式
- 所有其他欄位都會保留在結果中

---

## 📊 輸出結果

程式會產生一個Excel檔案（檔名_adduct_results.xlsx），包含**單一工作表**，並在原始數據上直接標記配對結果：

### Original_Data_Annotated（標記後的原始數據）

**新增的欄位：**
- **Adduct_Type**: 加合物類型（[M+H]+ 或其他加合物）
  - 如果有多個可能，用分號`;`分隔
  - 例如：`[M+Na]+; [M+K]+`
- **Description**: 配對說明
  - 例如：`[M+Na]+ of Base m/z 242.1191`
  - 多個配對：`[M+Na]+ of Base m/z 242.1000; [M+K]+ of Base m/z 220.0000`
- **Pair_mz**: 配對的Base m/z值
  - 單一配對：`242.1191`
  - 多個配對：`242.1000; 220.0000`（用分號分隔）
- **PPM_Error**: PPM誤差
  - 單一配對：`3.50`
  - 多個配對：`3.00; 5.00`（用分號分隔，對應Pair_mz的順序）

**視覺化標記：**
- ⚫ **黑色字體** = [M+H]+ 且未被配對（沒有找到對應的加合物）
- 🔴 **紅色字體** = [M+H]+ 且有被配對（是某個加合物的Base）
- 🟡 **黃色背景** = 非[M+H]+的加合物（如[M+Na]+, [M+K]+等）

**範例：**
假設有以下數據：
- m/z 150.0 → 未配對，顯示為**黑色字體**
- m/z 242.1191 (ID: 10) 配對到 264.1010 → 顯示為**紅色字體**（有配對的Base）
- m/z 264.1010 (ID: 15) 是 [M+Na]+ → 顯示為**黃色背景**

**多重配對範例：**
如果一個Pair可能來自多個Base：
- m/z 264.1 可能是：
  - m/z 242.1 的 [M+Na]+（PPM: 3.0）
  - m/z 220.0 的 [M+K]+（PPM: 5.0）
- 各欄位顯示：
  - `Adduct_Type`: `[M+Na]+; [M+K]+`
  - `Description`: `[M+Na]+ of Base m/z 242.1000; [M+K]+ of Base m/z 220.0000`
  - `Pair_mz`: `242.1000; 220.0000`
  - `PPM_Error`: `3.00; 5.00`
- 所有資訊按PPM誤差排序，最可能的配對列在前面

在Excel中你可以立即看出：
- 黑字 = 沒配對到加合物的[M+H]+
- 紅字 = 有配對到加合物的[M+H]+（Base）
- 黃底 = 這個訊號是加合物（可能來自多個Base）

**優點：**
- ✅ 一目了然，不需要對照多個工作表
- ✅ 保留所有原始欄位和數據
- ✅ 可以直接在原始數據上篩選和排序
- ✅ 輕鬆識別哪些訊號是加合物，可以選擇移除或保留

---

## ⚙️ 參數說明

### ppm_tolerance (PPM容許誤差)

**預設值**: 20 ppm

**計算方式**: 
- 使用** m/z 值**作為參考
- 容許範圍 (Da) = ppm × 較大的m/z值 / 1,000,000
- 例如：m/z 264.1010，20 ppm = 0.00528 Da 的容許範圍



**實際範例**:
```
設定 20 ppm:
- 對於 m/z 264.1010，容許範圍 = ±0.00528 Da
- 對於 m/z 500.0，容許範圍 = ±0.01000 Da
- 較大的 m/z 有較大的容許範圍（符合儀器特性）
```


### rt_tolerance (RT容許誤差)

**預設值**: 0.05 分鐘



## 💡 使用範例

### 範例1: 使用圖形介面（最簡單）

```bash
# 直接執行程式
python adduct_matcher.py

# 在彈出的視窗中:
# 1. 點擊「選擇檔案」
# 2. 選擇你的數據檔案
# 3. 設定參數（可選）
# 4. 點擊「開始比對」
```

### 範例2: 基本使用

```python
from adduct_matcher import AdductMatcher

# 建立比對器
matcher = AdductMatcher(ppm_tolerance=20.0)

# 處理檔案
results = matcher.process(
    file_path='LC_MS_data.xlsx',
    rt_tolerance=0.05
)

# 查看找到多少配對
print(f"找到 {len(results)} 個配對")
```

### 範例2: 基本使用

```python
from adduct_matcher import AdductMatcher

# 建立比對器
matcher = AdductMatcher(ppm_tolerance=20.0)

# 處理檔案
results = matcher.process(
    file_path='LC_MS_data.xlsx',
    rt_tolerance=0.05
)

# 查看找到多少配對
print(f"找到 {len(results)} 個配對")
```

### 範例3: 自訂輸出檔名

```python
matcher = AdductMatcher(ppm_tolerance=20.0)

results = matcher.process(
    file_path='data.xlsx',
    rt_tolerance=0.05,
    output_file='my_results.xlsx'
)
```

### 範例3: 高解析度質譜數據

```python
# 使用更嚴格的參數
matcher = AdductMatcher(ppm_tolerance=5.0)

results = matcher.process(
    file_path='HRMS_data.xlsx',
    rt_tolerance=0.02
)
```

### 範例3: 自訂輸出檔名

```python
matcher = AdductMatcher(ppm_tolerance=20.0)

results = matcher.process(
    file_path='data.xlsx',
    rt_tolerance=0.05,
    output_file='my_results.xlsx'
)
```

### 範例4: 高解析度質譜數據

```python
# 使用更嚴格的參數
matcher = AdductMatcher(ppm_tolerance=5.0)

results = matcher.process(
    file_path='HRMS_data.xlsx',
    rt_tolerance=0.02
)
```

### 範例5: 分步處理

```python
matcher = AdductMatcher(ppm_tolerance=20.0)

# 1. 載入數據
df = matcher.load_data('data.xlsx')

# 2. 執行配對
results = matcher.match_adducts(df, rt_tolerance=0.05)

# 3. 儲存結果
matcher.save_results(df, 'data.xlsx', results)
```

### 範例6: 批次處理

```python
import glob
from adduct_matcher import AdductMatcher

matcher = AdductMatcher(ppm_tolerance=20.0)

for file in glob.glob('*.xlsx'):
    print(f"\n處理: {file}")
    results = matcher.process(file, rt_tolerance=0.05)
    print(f"完成: 找到 {len(results)} 個配對")
```

---

## 🎓 實際應用案例

### 你的原始需求案例

**輸入數據**:
- RT = 10.99 分鐘
- m/z = 242.1191 (ID: 10)
- m/z = 264.1010 (ID: 15)

**程式處理**:
1. ✓ 識別兩者在相同RT（差異 < 0.05分鐘）
2. ✓ 計算質量差: 264.1010 - 242.1191 = 21.9819 Da
3. ✓ 比對加合物表，找到 [M+Na]+ 理論值 = 21.981944 Da
4. ✓ 計算PPM誤差 = 2.01 ppm（< 20 ppm ✓）
5. ✓ 標註: m/z 264.1010 是 ID 10 的 [M+Na]+ 加合物

**輸出結果**:
```
Base_Peak ID: 10
Base_m/z: 242.1191
Pair_Peak ID: 15
Pair_m/z: 264.1010
Pair_Adduct: [M+Na]+
PPM_Error: 2.01
(以及所有其他原始欄位)
```

---

## 🔧 常見問題排解

### Q1: 程式無法識別我的欄位

**可能原因**: 欄位名稱不在支援列表中

**解決方法**:
1. 查看上方"自動欄位識別"部分的支援名稱
2. 將欄位重命名為標準名稱（例如: RT, m/z, Intensity）
3. 或修改 `_find_column` 方法添加你的欄位名稱

### Q2: 沒有找到任何配對

**可能原因與解決方法**:

1. **RT容許範圍太小**
   ```python
   # 增加RT容許範圍
   results = matcher.process('data.xlsx', rt_tolerance=0.2)
   ```

2. **PPM容許範圍太嚴格**
   ```python
   # 增加PPM容許範圍
   matcher = AdductMatcher(ppm_tolerance=50.0)
   ```

3. **數據中沒有真實的加合物**
   - 檢查原始數據，確認是否有相似RT的訊號
   - 確認質量差是否接近表中的加合物

### Q3: 配對結果太多

**解決方法**:

1. **使用更嚴格的參數**
   ```python
   matcher = AdductMatcher(ppm_tolerance=5.0)
   results = matcher.process('data.xlsx', rt_tolerance=0.02)
   ```

2. **後處理篩選**
   ```python
   # 只保留PPM誤差 < 5 的結果
   high_quality = results[results['PPM_Error'] < 5.0]
   
   # 只看特定加合物
   na_only = results[results['Pair_Adduct'] == '[M+Na]+']
   ```

**注意**: 程式已自動實施策略A，當多個訊號競爭配對時，會自動保留RT最接近的配對。

### Q4: 如果有同分異構物怎麼辦？

如果你的樣品中可能有同分異構物（相同m/z但不同RT），程式會正確處理：
- 每個同分異構物會被視為獨立的[M+H]+
- 如果它們都有對應的加合物，都會被標記
- 程式使用m/z和RT雙重判斷，不會混淆不同RT的訊號

### Q4: 結果中看不到我的其他欄位

確認:
- 原始檔案中確實有這些欄位
- 檔案讀取成功（查看console輸出的"保留其他欄位"數量）
- 在結果檔案的 Adduct_Pairs 工作表中查找 Base_欄位名 和 Pair_欄位名

---

## 📈 支援的23種ESI加合物

| From | To | Delta (Da) |
|------|-----|-----------|
| [M+H]+ | [M+Li]+ | 6.00817839 |
| [M+H]+ | [M+NH4]+ | 17.02654909 |
| [M+H]+ | [M+Na]+ | 21.98194424 |
| [M+H]+ | [M+K]+ | 37.95588144 |
| [M+H]+ | [M+H3O]+ | 18.01056467 |
| [M+H]+ | [M+H2O+H]+ | 18.01056468 |
| [M+H]+ | [M+MeOH+H]+  | 32.02621475 |
| [M+H]+ | [M+EtOH+H]+ | 46.04186481 |
| [M+H]+ | [M+IPA+H]+ | 60.05751488 |
| [M+H]+ | [M+ACN+H]+ | 41.0265491 |
| [M+H]+ | [M+DMSO+H]+ | 78.01393599 |
| [M+H]+ | [M+H2O+Na]+ | 39.99250892 |
| [M+H]+ | [M+MeOH+Na]+ | 54.00815898 |
| [M+H]+ | [M+EtOH+Na]+ | 68.02380905 |
| [M+H]+ | [M+IPA+Na]+ | 82.03945911 |
| [M+H]+ | [M+ACN+Na]+ | 63.00849334 |
| [M+H]+ | [M+DMSO+Na]+ | 99.99588022 |
| [M+H]+ | [M+HCOOH+H]+ | 46.0054793 |
| [M+H]+ | [M+CH3COOH+H]+ | 60.02112937 |
| [M+H]+ | [M+H2CO3+H]+ | 62.00039392 |
| [M+H]+ | [M+H2SO4+H]+ | 97.96737972 |
| [M+H]+ | [M+Na+H]+ | 22.9892207 |
| [M+H]+ | [M+2Na-H]+ | 43.96388847 |

---

## 💻 程式架構

### 主要類別: AdductMatcher

**初始化**:
```python
AdductMatcher(ppm_tolerance=20.0)
```

**主要方法**:

1. `load_data(file_path)`: 載入數據並自動識別欄位
2. `match_adducts(df, rt_tolerance)`: 執行加合物比對
3. `save_results(df, original_file, results, output_file)`: 儲存結果
4. `process(file_path, rt_tolerance, output_file)`: 完整處理流程（推薦）

---

## 🎯 最佳實踐

1. **數據預處理**: 確保數據品質良好，移除明顯的雜訊峰
2. **參數優化**: 先用預設參數，根據結果調整
3. **結果驗證**: 檢查PPM誤差，通常好的配對在5 ppm以內
4. **保留原始數據**: 程式會自動保留所有欄位，方便後續分析
5. **批次分析**: 使用相同參數處理同一批樣品以保持一致性

---

## 📝 注意事項

1. 程式假設較小的m/z為 [M+H]+，較大的m/z為其他加合物
2. 只移除 m/z ≤ 0 或 Intensity ≤ 0 的數據，RT允許為0
3. 同一組配對只保留一次（根據 Base_m/z, Pair_m/z, Pair_Adduct去重）
4. 所有原始欄位都會被保留並加上 Base_ 或 Pair_ 前綴

---

## 📄 授權資訊

此程式由 Claude (Anthropic) 協助開發  
專為質譜數據分析設計
