#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ESI加合物比對程式
自動識別質譜數據中的加合物配對
"""

import pandas as pd
import numpy as np
from pathlib import Path
from typing import Optional, List, Dict, Tuple
import warnings
warnings.filterwarnings('ignore')


class Adduct_matcher:
    """ESI加合物比對器"""
    
    def __init__(self, ppm_tolerance: float = 20.0):
        """
        初始化加合物比對器
        
        Parameters:
        -----------
        ppm_tolerance : float
            質量誤差容許範圍 (ppm)，預設為 20 ppm
        """
        self.ppm_tolerance = ppm_tolerance / 1_000_000  # 轉換為比例
        self.adduct_table = self._create_adduct_table()
        
        # 欄位相關資訊
        self.rt_col = None
        self.mz_col = None
        self.intensity_col = None
        self.all_columns = None
        
    def _create_adduct_table(self) -> pd.DataFrame:
        """建立23個常見ESI加合物質量差異表"""
        data = {
            'From': ['[M+H]+'] * 23,
            'To': [
                '[M+Li]+', '[M+NH4]+', '[M+Na]+', '[M+K]+', '[M+H3O]+',
                '[M+H2O+H]+', '[M+MeOH+H]+', '[M+EtOH+H]+', '[M+IPA+H]+',
                '[M+ACN+H]+', '[M+DMSO+H]+', '[M+H2O+Na]+', '[M+MeOH+Na]+',
                '[M+EtOH+Na]+', '[M+IPA+Na]+', '[M+ACN+Na]+', '[M+DMSO+Na]+',
                '[M+HCOOH+H]+', '[M+CH3COOH+H]+', '[M+H2CO3+H]+',
                '[M+H2SO4+H]+', '[M+Na+H]+', '[M+2Na-H]+'
            ],
            'Delta_Da': [
                6.00817839, 17.02654909, 21.98194424, 37.95588144, 18.01056467,
                18.01056468, 32.02621475, 46.04186481, 60.05751488, 41.0265491,
                78.01393599, 39.99250892, 54.00815898, 68.02380905, 82.03945911,
                63.00849334, 99.99588022, 46.0054793, 60.02112937, 62.00039392,
                97.96737972, 22.9892207, 43.96388847
            ]
        }
        return pd.DataFrame(data)
    
    def load_data(self, file_path: str) -> pd.DataFrame:
        """
        載入數據並自動識別欄位 (支援 Excel, CSV, TSV)
        
        Parameters:
        -----------
        file_path : str
            檔案路徑
            
        Returns:
        --------
        pd.DataFrame
            包含所有欄位的數據框
        """
        file_path = str(file_path)
        
        # 根據副檔名選擇讀取方式
        if file_path.endswith('.csv'):
            df = pd.read_csv(file_path, encoding='utf-8-sig')
        elif file_path.endswith('.tsv') or file_path.endswith('.txt'):
            df = pd.read_csv(file_path, sep='\t', encoding='utf-8-sig')
        elif file_path.endswith(('.xlsx', '.xls', '.xlsm', '.xlsb')):
            df = pd.read_excel(file_path)
        else:
            raise ValueError(f"不支援的檔案格式。支援: .xlsx, .xls, .xlsm, .xlsb, .csv, .tsv, .txt")
        
        # 自動識別 RT, m/z, Intensity 欄位
        rt_col = self._find_column(df.columns, [
            'rt', 'retention time', 'retention_time', 'retentiontime',
            'rt (min)', 'rt(min)', 'retention time (min)', 'time',
            'r.t.', 'r.t', '保留時間', '滯留時間'
        ])
        mz_col = self._find_column(df.columns, [
            'm/z', 'mz', 'm_z', 'mass', 'mass to charge',
            'precursor ion m/z', 'precursor m/z', 'precursormz',
            '質量', '質荷比'
        ])
        intensity_col = self._find_column(df.columns, [
            'intensity', 'int', 'intens', 'abundance', 'abund', 'height', 'area',
            'precursor ion intensity', 'precursor intensity', 'precursorintensity',
            '強度', '峰高', '峰面積'
        ])
        
        if not rt_col or not mz_col or not intensity_col:
            missing = []
            if not rt_col: missing.append("RT (滯留時間)")
            if not mz_col: missing.append("m/z (質荷比)")
            if not intensity_col: missing.append("Intensity (強度)")
            
            available_cols = "\n可用的欄位: " + ", ".join(df.columns.tolist())
            raise ValueError(f"無法識別欄位: {', '.join(missing)}\n請確認標頭包含這些欄位名稱{available_cols}")
        
        # 標記主要欄位
        self.rt_col = rt_col
        self.mz_col = mz_col
        self.intensity_col = intensity_col
        self.all_columns = list(df.columns)
        
        print(f"✓ 成功讀取檔案: {Path(file_path).name}")
        print(f"  數據形狀: {df.shape[0]} 行 × {df.shape[1]} 欄")
        print(f"\n已識別欄位:")
        print(f"  RT (滯留時間):  {rt_col}")
        print(f"  m/z (質荷比):    {mz_col}")
        print(f"  Intensity (強度): {intensity_col}")
        print(f"  保留其他欄位: {len(self.all_columns) - 3} 個")
        
        # 移除無效數據 (只檢查 m/z 和 intensity > 0, RT 允許為 0)
        original_count = len(df)
        df = df[(df[mz_col] > 0) & (df[intensity_col] > 0)]
        df = df.dropna(subset=[rt_col, mz_col, intensity_col])
        
        removed_count = original_count - len(df)
        if removed_count > 0:
            print(f"\n⚠ 已移除 {removed_count} 行無效數據 (m/z≤0, Intensity≤0, 或包含缺失值)")
        
        print(f"✓ 有效數據: {len(df)} 行\n")
        
        return df.reset_index(drop=True)
    
    def _find_column(self, columns: List[str], possible_names: List[str]) -> Optional[str]:
        """
        尋找符合的欄位名稱
        
        Parameters:
        -----------
        columns : list
            所有欄位名稱
        possible_names : list
            可能的欄位名稱列表
            
        Returns:
        --------
        str or None
            找到的欄位名稱
        """
        for col in columns:
            col_lower = str(col).lower().strip()
            for name in possible_names:
                if name in col_lower:
                    return col
        return None
    
    def match_adducts(self, df: pd.DataFrame, rt_tolerance: float = 0.05) -> pd.DataFrame:
        """
        比對加合物，保留所有原始欄位
        當多個base競爭同一個pair時，保留所有配對（列出所有可能）
        
        Parameters:
        -----------
        df : pd.DataFrame
            包含質譜數據的DataFrame
        rt_tolerance : float
            滯留時間容許誤差（分鐘），預設為0.05分鐘
            
        Returns:
        --------
        pd.DataFrame
            包含加合物比對結果的DataFrame，保留所有原始欄位
        """
        print("開始進行加合物比對...")
        
        # 確保數據按RT排序
        df = df.sort_values(by=self.rt_col).reset_index(drop=True)
        
        # 建立結果列表
        results = []
        
        # 對每個數據點進行比對
        for i, row in df.iterrows():
            current_rt = row[self.rt_col]
            current_mz = row[self.mz_col]
            current_intensity = row[self.intensity_col]
            
            # 找出相同或接近滯留時間的其他訊號
            rt_mask = abs(df[self.rt_col] - current_rt) <= rt_tolerance
            nearby_peaks = df[rt_mask & (df.index != i)]
            
            # 與附近的peak比對
            for _, nearby_peak in nearby_peaks.iterrows():
                nearby_rt = nearby_peak[self.rt_col]
                nearby_mz = nearby_peak[self.mz_col]
                nearby_intensity = nearby_peak[self.intensity_col]
                
                # 計算RT差異（用於後續篩選）
                rt_diff = abs(current_rt - nearby_rt)
                
                # 計算質量差 (大減小)
                if nearby_mz > current_mz:
                    mz_diff = nearby_mz - current_mz
                    base_mz = current_mz
                    base_rt = current_rt
                    base_row = row
                    pair_mz = nearby_mz
                    pair_rt = nearby_rt
                    pair_row = nearby_peak
                else:
                    mz_diff = current_mz - nearby_mz
                    base_mz = nearby_mz
                    base_rt = nearby_rt
                    base_row = nearby_peak
                    pair_mz = current_mz
                    pair_rt = current_rt
                    pair_row = row
                
                # 與加合物表比對
                for _, adduct in self.adduct_table.iterrows():
                    theoretical_delta = adduct['Delta_Da']
                    
                    # 計算ppm容許範圍 (使用較大的m/z作為參考，與VBA相同)
                    reference_mz = max(current_mz, nearby_mz)  # 使用較大的m/z
                    ppm_tolerance_da = self.ppm_tolerance * reference_mz  # 轉換為Da
                    
                    # 檢查質量差是否在容許範圍內
                    if abs(mz_diff - theoretical_delta) <= ppm_tolerance_da:
                        # 計算實際的PPM誤差（用於報告）
                        ppm_error_value = abs(mz_diff - theoretical_delta) / reference_mz * 1_000_000
                        
                        # 建立結果字典，包含所有原始欄位
                        result = {}
                        
                        # 添加Base化合物的所有欄位（加上Base_前綴）
                        for col in self.all_columns:
                            result[f'Base_{col}'] = base_row[col]
                        
                        # 添加Pair化合物的所有欄位（加上Pair_前綴）
                        for col in self.all_columns:
                            result[f'Pair_{col}'] = pair_row[col]
                        
                        # 添加比對資訊
                        result['Base_Adduct'] = '[M+H]+'
                        result['Pair_Adduct'] = adduct['To']
                        result['Theoretical_Delta_Da'] = theoretical_delta
                        result['Observed_Delta_Da'] = mz_diff
                        result['PPM_Error'] = round(ppm_error_value, 2)
                        result['RT_Diff'] = rt_diff  # 添加RT差異用於篩選
                        result['Reference_mz'] = reference_mz
                        result['Annotation'] = f"{adduct['To']} of Base (m/z {base_mz:.4f})"
                        
                        results.append(result)
        
        if results:
            results_df = pd.DataFrame(results)
            
            # 移除重複的配對（同一組Base-Pair-Adduct只保留一次）
            base_mz_col = f'Base_{self.mz_col}'
            base_rt_col = f'Base_{self.rt_col}'
            pair_mz_col = f'Pair_{self.mz_col}'
            pair_rt_col = f'Pair_{self.rt_col}'
            
            # 先去除完全相同的配對
            results_df = results_df.drop_duplicates(
                subset=[base_mz_col, base_rt_col, pair_mz_col, pair_rt_col, 'Pair_Adduct'], 
                keep='first'
            ).reset_index(drop=True)
            
            # 方案2: 保留所有配對，但在同一個Base有多個Pair時，選RT最近的
            # 當一個Base有多個Pair候選時，保留RT最接近的
            results_df = results_df.sort_values('RT_Diff')
            results_df = results_df.drop_duplicates(
                subset=[base_mz_col, base_rt_col, 'Pair_Adduct'],
                keep='first'
            ).reset_index(drop=True)
            
            # 當一個Pair有多個Base時，保留所有配對（不去重）
            # 這樣可以列出所有可能性
            
            # 移除RT_Diff欄位（用戶不需要看到）
            results_df = results_df.drop(columns=['RT_Diff'])
            
            print(f"✓ 找到 {len(results_df)} 個加合物配對")
            
            # 顯示加合物類型統計
            adduct_counts = results_df['Pair_Adduct'].value_counts()
            print(f"\n加合物類型分布:")
            for adduct, count in adduct_counts.head(5).items():
                print(f"  {adduct}: {count} 個")
            
            return results_df
        else:
            print("⚠ 未找到符合的加合物配對")
            return pd.DataFrame()
    
    def save_results(self, df: pd.DataFrame, original_file: str, 
                    results: pd.DataFrame, output_file: Optional[str] = None):
        """
        儲存結果到Excel檔案，在原始數據上標記配對結果
        
        Parameters:
        -----------
        df : pd.DataFrame
            原始數據
        original_file : str
            原始檔案路徑
        results : pd.DataFrame
            配對結果
        output_file : str, optional
            輸出檔案名稱，若未指定則自動生成
        """
        if output_file is None:
            # 自動生成輸出檔案名稱
            file_path = Path(original_file)
            output_file = str(file_path.parent / f"{file_path.stem}_adduct_results.xlsx")
        
        print(f"\n正在儲存結果到: {Path(output_file).name}")
        
        # 準備標記資訊
        df_marked = df.copy()
        
        # 新增識別欄位
        df_marked['Adduct_Type'] = '[M+H]+'  # 預設都是[M+H]+
        df_marked['Pair_mz'] = ''            # 配對的m/z值
        df_marked['PPM_Error'] = ''          # PPM誤差
        df_marked['Description'] = ''     # 配對說明
        df_marked['Is_Matched_Base'] = False  # 是否為有配對的Base
        
        if not results.empty:
            # 建立m/z到索引的映射
            mz_col = self.mz_col
            rt_col = self.rt_col
            base_mz_col = f'Base_{mz_col}'
            pair_mz_col = f'Pair_{mz_col}'
            base_rt_col = f'Base_{rt_col}'
            pair_rt_col = f'Pair_{rt_col}'
            
            # 建立Pair到多個Base的映射字典
            pair_to_bases = {}  # key: (pair_mz, pair_rt), value: list of (base_mz, adduct, ppm)
            
            for _, result_row in results.iterrows():
                base_mz = result_row[base_mz_col]
                base_rt = result_row[base_rt_col]
                pair_mz = result_row[pair_mz_col]
                pair_rt = result_row[pair_rt_col]
                adduct = result_row['Pair_Adduct']
                ppm = result_row['PPM_Error']
                
                pair_key = (pair_mz, pair_rt)
                if pair_key not in pair_to_bases:
                    pair_to_bases[pair_key] = []
                
                pair_to_bases[pair_key].append({
                    'base_mz': base_mz,
                    'adduct': adduct,
                    'ppm': ppm
                })
            
            # 標記Base化合物（有配對的[M+H]+）
            for _, result_row in results.iterrows():
                base_mz = result_row[base_mz_col]
                base_rt = result_row[base_rt_col]
                
                # 同時檢查 m/z 和 RT
                mask = (abs(df_marked[mz_col] - base_mz) < 0.0001) & \
                       (abs(df_marked[rt_col] - base_rt) < 0.0001)
                if mask.any():
                    df_marked.loc[mask, 'Is_Matched_Base'] = True
            
            # 標記Pair化合物（標記為對應的加合物類型，可能有多個）
            for pair_key, base_list in pair_to_bases.items():
                pair_mz, pair_rt = pair_key
                
                # 找到對應的Pair行
                mask = (abs(df_marked[mz_col] - pair_mz) < 0.0001) & \
                       (abs(df_marked[rt_col] - pair_rt) < 0.0001)
                
                if mask.any():
                    # 如果有多個Base，按PPM排序（最小的最可能）
                    base_list_sorted = sorted(base_list, key=lambda x: x['ppm'])
                    
                    # 合併所有加合物類型
                    adduct_types = '; '.join([b['adduct'] for b in base_list_sorted])
                    
                    # 合併所有配對說明（不包含PPM資訊）
                    descriptions = '; '.join([
                        f"{b['adduct']} of Base m/z {b['base_mz']:.4f}"
                        for b in base_list_sorted
                    ])
                    
                    # 合併所有可能的 Base m/z（用分號分隔）
                    pair_mz_list = '; '.join([f"{b['base_mz']:.4f}" for b in base_list_sorted])
                    
                    # 合併所有PPM誤差（用分號分隔）
                    ppm_list = '; '.join([f"{b['ppm']:.2f}" for b in base_list_sorted])
                    
                    df_marked.loc[mask, 'Adduct_Type'] = adduct_types
                    df_marked.loc[mask, 'Pair_mz'] = pair_mz_list  # 列出所有Base m/z
                    df_marked.loc[mask, 'PPM_Error'] = ppm_list  # 列出所有PPM誤差
                    df_marked.loc[mask, 'Description'] = descriptions
        
        # 寫入Excel並設定格式
        with pd.ExcelWriter(output_file, engine='openpyxl') as writer:
            # 寫入標記後的原始數據（不包含Is_Matched_Base欄位）
            df_output = df_marked.drop(columns=['Is_Matched_Base'])
            df_output.to_excel(writer, sheet_name='Original_Data_Annotated', index=False)
            
            # 取得workbook和worksheet以設定格式
            from openpyxl.styles import PatternFill, Font
            workbook = writer.book
            worksheet = writer.sheets['Original_Data_Annotated']
            
            if not results.empty:
                # 找到Adduct_Type欄位的索引
                adduct_type_col_idx = list(df_output.columns).index('Adduct_Type') + 1
                
                # 設定顏色
                yellow_fill = PatternFill(start_color='FFFF00', end_color='FFFF00', fill_type='solid')
                red_font = Font(color='FF0000', bold=True)
                
                # 遍歷每一行設定格式
                for row_idx in range(2, len(df_marked) + 2):  # 從第2行開始（第1行是標題）
                    adduct_type_value = df_marked.iloc[row_idx-2]['Adduct_Type']
                    is_matched_base = df_marked.iloc[row_idx-2]['Is_Matched_Base']
                    
                    if adduct_type_value == '[M+H]+' and is_matched_base:
                        # 有配對的Base化合物 - 紅色字體
                        for col_idx in range(1, len(df_output.columns) + 1):
                            cell = worksheet.cell(row=row_idx, column=col_idx)
                            cell.font = red_font
                    elif adduct_type_value and adduct_type_value != '[M+H]+':
                        # Pair化合物 - 黃色背景
                        for col_idx in range(1, len(df_output.columns) + 1):
                            cell = worksheet.cell(row=row_idx, column=col_idx)
                            cell.fill = yellow_fill
                    # 未配對的[M+H]+ 保持黑色（不做任何處理）
                
                # 格式化 Intensity 欄位為科學記號
                intensity_col_idx = list(df_marked.columns).index(self.intensity_col) + 1
                for row in range(2, len(df_marked) + 2):
                    cell = worksheet.cell(row=row, column=intensity_col_idx)
                    cell.number_format = '0.00E+00'
                
                # 如果有PPM_Error欄位，也格式化
                if 'PPM_Error' in df_marked.columns:
                    ppm_col_idx = list(df_marked.columns).index('PPM_Error') + 1
                    for row in range(2, len(df_marked) + 2):
                        cell = worksheet.cell(row=row, column=ppm_col_idx)
                        if cell.value and cell.value != '':
                            cell.number_format = '0.00'
                
                # 調整欄寬
                for column in worksheet.columns:
                    max_length = 0
                    column_letter = column[0].column_letter
                    for cell in column:
                        try:
                            if len(str(cell.value)) > max_length:
                                max_length = len(str(cell.value))
                        except:
                            pass
                    adjusted_width = min(max_length + 2, 50)
                    worksheet.column_dimensions[column_letter].width = adjusted_width
            else:
                # 沒有配對結果時的簡單格式化
                df_marked.to_excel(writer, sheet_name='Original_Data_Annotated', index=False)
        
        print(f"✓ 結果已成功儲存!")
        print(f"  Original_Data_Annotated - 標記後的原始數據")
        if not results.empty:
            base_count = df_marked['Is_Matched_Base'].sum()
            pair_count = (df_marked['Adduct_Type'] != '[M+H]+').sum()
            unmatched_count = len(df_marked) - base_count - pair_count
            
            print(f"    • 黑色字體 = 未配對的 [M+H]+ ({unmatched_count} 個)")
            print(f"    • 紅色字體 = 有配對的 Base化合物 ([M+H]+) ({base_count} 個)")
            print(f"    • 黃色背景 = 配對的加合物 ({pair_count} 個)")
            print(f"    • 新增欄位: Adduct_Type, Description, Pair_mz, PPM_Error")
        else:
            print(f"    • 未找到配對結果")
    
    def process(self, file_path: str, rt_tolerance: float = 0.05, 
                output_file: Optional[str] = None) -> pd.DataFrame:
        """
        完整處理流程
        
        Parameters:
        -----------
        file_path : str
            輸入檔案路徑
        rt_tolerance : float
            RT容許誤差（分鐘）
        output_file : str, optional
            輸出檔案名稱
            
        Returns:
        --------
        pd.DataFrame
            配對結果
        """
        print("="*70)
        print("ESI加合物比對程式")
        print("="*70 + "\n")
        
        # 載入數據
        df = self.load_data(file_path)
        
        # 顯示數據預覽
        print("數據預覽 (前5行):")
        print(df.head().to_string(index=False))
        print()
        
        # 執行配對
        results = self.match_adducts(df, rt_tolerance=rt_tolerance)
        
        # 儲存結果
        if not results.empty:
            self.save_results(df, file_path, results, output_file)
            
            print("\n" + "="*70)
            print("✓ 分析完成!")
            print("="*70)
            
            # 顯示關鍵統計
            print(f"\n找到 {len(results)} 個加合物配對")
            print(f"平均PPM誤差: {results['PPM_Error'].mean():.2f}")
            print(f"PPM誤差範圍: {results['PPM_Error'].min():.2f} - {results['PPM_Error'].max():.2f}")
        else:
            print("\n" + "="*70)
            print("⚠ 未找到符合的加合物配對")
            print("="*70)
            print("\n建議:")
            print("  1. 增加RT容許誤差 (例如: 0.1 或 0.2)")
            print("  2. 增加PPM容許誤差 (例如: 30 或 50)")
            print("  3. 檢查數據品質")
        
        return results


class Adduct_matcherGUI:
    """圖形化介面"""
    
    def __init__(self, root):
        import tkinter as tk
        from tkinter import filedialog, messagebox
        
        self.tk = tk
        self.filedialog = filedialog
        self.messagebox = messagebox
        
        self.root = root
        self.root.title("ESI加合物比對程式")
        self.root.geometry("600x500")
        
        self.input_file = None
        
        self.create_widgets()
    
    def create_widgets(self):
        tk = self.tk
        
        # 標題
        title_label = tk.Label(
            self.root, 
            text="ESI加合物比對程式",
            font=("Arial", 16, "bold")
        )
        title_label.pack(pady=10)
        
        # 檔案選擇框架
        file_frame = tk.LabelFrame(self.root, text="檔案選擇", padx=10, pady=10)
        file_frame.pack(padx=20, pady=10, fill="x")
        
        self.file_label = tk.Label(file_frame, text="未選擇檔案", fg="gray")
        self.file_label.pack(side="left", padx=5)
        
        tk.Button(
            file_frame, 
            text="選擇檔案 (Excel/CSV/TSV)", 
            command=self.select_file
        ).pack(side="right", padx=5)
        
        # 參數設定框架
        param_frame = tk.LabelFrame(self.root, text="參數設定", padx=10, pady=10)
        param_frame.pack(padx=20, pady=10, fill="x")
        
        # PPM 容差
        tk.Label(param_frame, text="PPM 容差:").grid(row=0, column=0, sticky="w", pady=5)
        self.ppm_tolerance_var = tk.StringVar(value="20")
        tk.Entry(param_frame, textvariable=self.ppm_tolerance_var, width=15).grid(row=0, column=1, pady=5)
        
        # RT 容差
        tk.Label(param_frame, text="RT 容差 (分鐘):").grid(row=1, column=0, sticky="w", pady=5)
        self.rt_tolerance_var = tk.StringVar(value="0.05")
        tk.Entry(param_frame, textvariable=self.rt_tolerance_var, width=15).grid(row=1, column=1, pady=5)
        
        # 執行按鈕
        tk.Button(
            self.root, 
            text="開始比對", 
            command=self.process_data,
            bg="#4CAF50",
            fg="white",
            font=("Arial", 12, "bold"),
            padx=20,
            pady=10
        ).pack(pady=20)
        
        # 狀態顯示
        self.status_text = tk.Text(self.root, height=12, width=70, state="disabled")
        self.status_text.pack(padx=20, pady=10)
    
    def select_file(self):
        """選擇輸入檔案"""
        file_path = self.filedialog.askopenfilename(
            title="選擇質譜數據檔案",
            filetypes=[
                ("所有支援格式", "*.xlsx *.xls *.xlsm *.xlsb *.csv *.tsv *.txt"),
                ("Excel files", "*.xlsx *.xls *.xlsm *.xlsb"),
                ("CSV files", "*.csv"),
                ("TSV files", "*.tsv *.txt"),
                ("All files", "*.*")
            ]
        )
        
        if file_path:
            self.input_file = file_path
            self.file_label.config(text=Path(file_path).name, fg="black")
    
    def update_status(self, message):
        """更新狀態顯示"""
        self.status_text.config(state="normal")
        self.status_text.insert("end", message + "\n")
        self.status_text.see("end")
        self.status_text.config(state="disabled")
        self.root.update()
    
    def process_data(self):
        """處理數據"""
        if not self.input_file:
            self.messagebox.showerror("錯誤", "請先選擇輸入檔案!")
            return
        
        try:
            # 清空狀態
            self.status_text.config(state="normal")
            self.status_text.delete(1.0, "end")
            self.status_text.config(state="disabled")
            
            # 讀取參數
            ppm_tol = float(self.ppm_tolerance_var.get())
            rt_tol = float(self.rt_tolerance_var.get())
            
            self.update_status("="*60)
            self.update_status("開始處理...")
            self.update_status("="*60)
            
            # 建立比對器
            matcher = Adduct_matcher(ppm_tolerance=ppm_tol)
            
            # 載入數據
            self.update_status("\n讀取數據中...")
            df = matcher.load_data(self.input_file)
            
            # 顯示識別的欄位
            self.update_status(f"\n已識別欄位:")
            self.update_status(f"  RT: {matcher.rt_col}")
            self.update_status(f"  m/z: {matcher.mz_col}")
            self.update_status(f"  Intensity: {matcher.intensity_col}")
            self.update_status(f"  保留其他欄位: {len(matcher.all_columns) - 3} 個")
            
            # 執行比對
            self.update_status("\n執行加合物比對...")
            results = matcher.match_adducts(df, rt_tolerance=rt_tol)
            
            if not results.empty:
                # 顯示加合物類型統計
                adduct_counts = results['Pair_Adduct'].value_counts()
                self.update_status(f"\n加合物類型分布:")
                for adduct, count in adduct_counts.head(5).items():
                    self.update_status(f"  {adduct}: {count} 個")
                
                # 生成輸出檔名
                input_path = Path(self.input_file)
                output_path = input_path.parent / f"{input_path.stem}_adduct_results.xlsx"
                
                # 儲存結果
                self.update_status("\n儲存結果中...")
                matcher.save_results(df, self.input_file, results, str(output_path))
                
                # 顯示統計
                self.update_status("\n" + "="*60)
                self.update_status("✓ 處理完成!")
                self.update_status("="*60)
                self.update_status(f"\n找到 {len(results)} 個加合物配對")
                self.update_status(f"平均PPM誤差: {results['PPM_Error'].mean():.2f}")
                self.update_status(f"PPM誤差範圍: {results['PPM_Error'].min():.2f} - {results['PPM_Error'].max():.2f}")
                self.update_status(f"\n結果已儲存至:\n{output_path}")
                
                self.messagebox.showinfo("完成", f"處理完成!\n\n找到 {len(results)} 個加合物配對\n\n結果已儲存至:\n{output_path}")
            else:
                self.update_status("\n" + "="*60)
                self.update_status("⚠ 未找到符合的加合物配對")
                self.update_status("="*60)
                self.update_status("\n建議:")
                self.update_status("  1. 增加RT容許誤差 (例如: 0.1 或 0.2)")
                self.update_status("  2. 增加PPM容許誤差 (例如: 30 或 50)")
                self.update_status("  3. 檢查數據品質")
                
                self.messagebox.showwarning("提醒", "未找到符合的加合物配對\n\n請嘗試調整參數或檢查數據")
            
        except Exception as e:
            self.messagebox.showerror("錯誤", f"處理時發生錯誤:\n{str(e)}")
            self.update_status(f"\n錯誤: {str(e)}")


def main():
    """主程式"""
    import tkinter as tk
    from tkinter import filedialog, messagebox
    
    root = tk.Tk()
    app = Adduct_matcherGUI(root)
    root.mainloop()


if __name__ == "__main__":
    main()
