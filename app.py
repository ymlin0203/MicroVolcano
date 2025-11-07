import streamlit as st
import pandas as pd
import numpy as np
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
import plotly.express as px

# ---------------------------
# Page & Style
# ---------------------------
st.set_page_config(page_title="🧬 Volcano Analyzer V10", layout="wide")
st.markdown("""
    <style>
    .main { background-color: #f8f9fa; }
    h1, h2, h3 { color: #2c3e50; font-family: 'Arial'; }
    .stButton>button {
        background-color: #4CAF50; color: white; border-radius: 8px;
        padding: 8px 16px; font-size: 16px;
    }
    </style>
""", unsafe_allow_html=True)

st.title("🧬 Volcano Analyzer V10 (2-Group Compare, FDR/p-value Toggle, Top5, Boxplot, Core)")

# ---------------------------
# Tabs（含 Core microbiota）
# ---------------------------
tab1, tab2, tab3, tab4, tab5, tab_core = st.tabs([
    "📂 上傳資料", "🌋 火山圖", "📊 結果表", "🏆 Top 5", "🧰 偵錯/Boxplot", "🌐 Core microbiota"
])

with tab1:
    st.subheader("📂 上傳資料")
    abundance_file = st.file_uploader("上傳 Genus abundance 表 (.tsv)", type=["tsv"])
    metadata_file = st.file_uploader("上傳 metadata 表 (.xlsx)", type=["xlsx"])

    # 分析參數
    log2fc_cutoff = st.slider("Log2 Fold Change 閾值", 0.5, 3.0, 1.0, 0.1)
    sig_mode = st.radio("顯著性判斷依據", ["FDR（Benjamini–Hochberg）", "p-value（不做 FDR）"], index=0)
    colA, colB = st.columns(2)
    with colA:
        if sig_mode.startswith("FDR"):
            fdr_cutoff = st.slider("FDR 閾值", 0.0, 0.2, 0.05, 0.005)
        else:
            p_cutoff = st.slider("p-value 閾值", 1e-4, 0.1, 0.05, 1e-4, format="%.4f")
    with colB:
        p_line = st.checkbox("在火山圖上顯示 p=0.05 水平線", value=True)

# ---------------------------
# Main logic
# ---------------------------
if abundance_file and metadata_file:
    # 讀 abundance（tsv），第一欄為 Genus，其餘為 Sample 列
    genus_df = pd.read_csv(abundance_file, sep="\t")
    genus_df = genus_df.rename(columns={genus_df.columns[0]: "Genus"})
    # 強制數值化
    for c in genus_df.columns[1:]:
        genus_df[c] = pd.to_numeric(genus_df[c], errors="coerce")

    # 讀 metadata
    meta_df = pd.read_excel(metadata_file, sheet_name=0)

    with tab1:
        # 選 SampleID 欄（避免抓錯第一欄）
        sample_col_in_meta = st.selectbox(
            "🔖 選擇 metadata 的 SampleID 欄位",
            list(meta_df.columns),
            index=(list(meta_df.columns).index("Sample") if "Sample" in meta_df.columns else 0)
        )
        # 產生/統一 SampleID 欄
        meta_df = meta_df.rename(columns={sample_col_in_meta: "SampleID"})
        meta_df["SampleID"] = meta_df["SampleID"].astype(str)

        # 分組欄位（排除 SampleID）
        candidate_cols = [c for c in meta_df.columns if c != "SampleID"]
        if not candidate_cols:
            st.error("metadata 檔沒有可用分組欄位。請至少提供 1 個分組欄。")
            st.stop()
        label_col = st.selectbox("🧭 分組欄位", candidate_cols, index=0)
        meta_df[label_col] = meta_df[label_col].astype(str)

        group_df = meta_df[["SampleID", label_col]].dropna()
        unique_vals = sorted(group_df[label_col].unique())
        if len(unique_vals) < 2:
            st.error("分組欄位只有一個值，無法做兩組比較。")
            st.stop()

        group_a = st.selectbox("🅰 Group A", unique_vals, index=0)
        group_b = st.selectbox("🅱 Group B", unique_vals, index=1 if len(unique_vals) > 1 else 0)

    if group_a == group_b:
        st.warning("⚠️ 請選擇兩個不同的 Group。")
        st.stop()

    # 僅保留同時存在於 abundance 的樣本
    sample_ids = [sid for sid in group_df["SampleID"] if sid in genus_df.columns]
    if len(sample_ids) == 0:
        st.error("❗ 找不到任何與 abundance 檔案相符的 SampleID，請確認兩檔案的樣本名稱一致。")
        st.stop()

    filtered_genus_df = genus_df[["Genus"] + sample_ids]

    # sample x genus 的矩陣
    data_matrix = filtered_genus_df.set_index("Genus").T
    # 若為 count，轉相對豐度
    needs_normalization = data_matrix.sum(axis=1).max() > 1.1
    if needs_normalization:
        st.info("⚙ 偵測到 count 資料，已自動轉為相對豐度。")
        data_matrix = data_matrix.div(data_matrix.sum(axis=1), axis=0)

    data_matrix.index.name = "SampleID"
    merged = data_matrix.merge(group_df, left_index=True, right_on="SampleID")
    merged = merged[merged[label_col].isin([group_a, group_b])]

    if merged.empty:
        st.error("❗ 合併後資料為空，請檢查分組欄位或 SampleID 是否一致。")
        st.stop()

    # ---------------------------
    # Statistics
    # ---------------------------
    results = []
    for genus in data_matrix.columns:
        g1 = merged.loc[merged[label_col] == group_a, genus].astype(float)
        g2 = merged.loc[merged[label_col] == group_b, genus].astype(float)

        # 資料不足跳過
        if g1.notna().sum() == 0 or g2.notna().sum() == 0:
            continue
        # 兩組皆零變異就跳過
        if (g1.std(ddof=0) == 0) and (g2.std(ddof=0) == 0):
            continue

        stat, p = mannwhitneyu(g1, g2, alternative="two-sided")
        # 注意：方向 A/B 會影響解讀，這裡維持原本 A/B
        log2_fc = np.log2((g1.mean() + 1e-6) / (g2.mean() + 1e-6))

        results.append({
            "Genus": genus,
            "p-value": float(p),
            "Group_A_Mean": float(g1.mean()),
            "Group_B_Mean": float(g2.mean()),
            "Log2_FC": float(log2_fc)
        })

    result_df = pd.DataFrame(results)
    if result_df.empty:
        st.warning("⚠️ 沒有可用的變數（可能兩組都為零或皆無變異）。")
        st.stop()

    # 統一衍生欄位
    result_df["-log10(p)"] = -np.log10(result_df["p-value"])

    # FDR 計算（即使之後不使用，表格也可顯示）
    _, fdr_p, _, _ = multipletests(result_df["p-value"], method='fdr_bh')
    result_df["FDR_p"] = fdr_p

    # 依據使用者選擇定義 Category
    eps = np.finfo(float).eps
    if sig_mode.startswith("FDR"):
        def _call(row):
            if (row["FDR_p"] <= fdr_cutoff + eps) and (row["Log2_FC"] > log2fc_cutoff):
                return "Upregulated"
            if (row["FDR_p"] <= fdr_cutoff + eps) and (row["Log2_FC"] < -log2fc_cutoff):
                return "Downregulated"
            return "Non-significant"
        caption = f"門檻：|Log2FC| > {log2fc_cutoff}，FDR ≤ {fdr_cutoff}"
        top_sort_key = "FDR_p"
    else:
        def _call(row):
            if (row["p-value"] <= p_cutoff + eps) and (row["Log2_FC"] > log2fc_cutoff):
                return "Upregulated"
            if (row["p-value"] <= p_cutoff + eps) and (row["Log2_FC"] < -log2fc_cutoff):
                return "Downregulated"
            return "Non-significant"
        caption = f"門檻：|Log2FC| > {log2fc_cutoff}，p ≤ {p_cutoff}"
        top_sort_key = "p-value"

    result_df["Category"] = result_df.apply(_call, axis=1)

    up_count = (result_df["Category"] == "Upregulated").sum()
    down_count = (result_df["Category"] == "Downregulated").sum()

    # ---------------------------
    # Volcano Plot
    # ---------------------------
    with tab2:
        st.subheader("🌋 Volcano Plot")
        st.caption(caption)
        fig = px.scatter(
            result_df,
            x="Log2_FC",
            y="-log10(p)",
            color="Category",
            hover_data=["Genus", "p-value", "FDR_p", "Group_A_Mean", "Group_B_Mean"],
            category_orders={"Category": ["Non-significant", "Upregulated", "Downregulated"]},
            color_discrete_map={
                "Non-significant": "gray",
                "Upregulated": "red",
                "Downregulated": "blue",
            }
        )
        if p_line:
            fig.add_hline(y=-np.log10(0.05), line_dash="dash")
        fig.add_vline(x=log2fc_cutoff, line_dash="dash")
        fig.add_vline(x=-log2fc_cutoff, line_dash="dash")
        fig.update_layout(legend_title_text="")  # 移除 legend 標題
        st.plotly_chart(fig, use_container_width=True)

        st.markdown(f"**Upregulated: {up_count}**　｜　**Downregulated: {down_count}**")

    # ---------------------------
    # 結果表 + 下載
    # ---------------------------
    with tab3:
        st.subheader("📊 結果表（含 p 與 FDR）")
        st.dataframe(
            result_df.sort_values(top_sort_key).style.apply(
                lambda s: ['background-color: #ffe6e6' if (s.name=='FDR_p' and v <= (fdr_cutoff if sig_mode.startswith('FDR') else 0)) else '' for v in s],
                subset=['FDR_p']
            ),
            use_container_width=True
        )
        st.download_button(
            "📥 下載結果 CSV",
            data=result_df.to_csv(index=False).encode(),
            file_name="2group_comparison.csv"
        )

    # ---------------------------
    # Top 5（依選定顯著性依據排序）
    # ---------------------------
    with tab4:
        st.subheader("🏆 Top 5")
        if sig_mode.startswith("FDR"):
            mask = (result_df["FDR_p"] <= fdr_cutoff + eps) & (result_df["Log2_FC"].abs() > log2fc_cutoff)
        else:
            mask = (result_df["p-value"] <= p_cutoff + eps) & (result_df["Log2_FC"].abs() > log2fc_cutoff)
        top5 = result_df[mask].sort_values(top_sort_key).head(5).copy()
        top5["Genus_name"] = top5["Genus"].apply(lambda x: x.split("__")[-1] if "__" in x else x)
        st.dataframe(top5[["Genus_name", "Log2_FC", "p-value", "FDR_p", "Category"]], use_container_width=True)
        st.download_button(
            "📥 下載 Top5 Summary",
            data=top5.to_csv(index=False).encode(),
            file_name="top5_summary.csv",
            disabled=top5.empty
        )

    # ---------------------------
    # 偵錯資訊 & Boxplot
    # ---------------------------
    with tab5:
        with st.expander("🔎 偵錯資訊（看不到圖時展開）"):
            st.write("group_df 頭：", group_df.head())
            st.write("data_matrix 形狀：", data_matrix.shape)
            st.write("merged 形狀：", merged.shape)
            st.write("result_df 形狀：", result_df.shape)
            st.write(result_df.head())

        st.subheader("📦 單一菌屬 Boxplot（選擇一個 Genus）")
        genus_choices = sorted(result_df["Genus"].unique().tolist())
        pick = st.selectbox("選擇 Genus", genus_choices)
        if pick:
            sub_box = merged[["SampleID", label_col, pick]].rename(columns={pick: "Abundance"})
            sub_box = sub_box.dropna(subset=["Abundance"])
            if sub_box.empty:
                st.info("此 Genus 在兩組皆無有效數據。")
            else:
                bfig = px.box(sub_box, x=label_col, y="Abundance", points="all", title=f"Boxplot: {pick}")
                st.plotly_chart(bfig, use_container_width=True)

    # ---------------------------
    # 🌐 Core microbiota 分頁
    # ---------------------------
    with tab_core:
        st.subheader("🌐 Core microbiota（核心菌群）")

        # 參數
        colA, colB, colC = st.columns(3)
        with colA:
            presence_thresh = st.number_input(
                "Presence 門檻（算“出現”的最低相對豐度）", min_value=0.0, max_value=0.05, value=0.001, step=0.0005,
                help="例如 0.001 = 0.1% 相對豐度以上才算有出現"
            )
        with colB:
            prev_cut = st.slider(
                "Prevalence 門檻（%樣本需出現）", min_value=0, max_value=100, value=80, step=5,
                help="例如 80 表示 ≥80% 樣本出現才算核心"
            )
        with colC:
            pooled_mean_min = st.number_input(
                "最小平均豐度（pooled mean）", min_value=0.0, max_value=0.05, value=0.01, step=0.001,
                help="例如 0.01 = 平均相對豐度 ≥ 1%"
            )

        mode = st.radio(
            "核心判定模式",
            ["Overall（整體）", "Per-group 交集（跨組皆達標）"],
            help="Overall：把所選組別的樣本合併計算一次 Prevalence；Per-group：各組都達標，再取交集"
        )

        # 可選納入哪些組別（預設全選）
        groups_all = sorted(merged[label_col].unique().tolist())
        groups_pick = st.multiselect(
            "納入的組別", groups_all, default=groups_all,
            help="只會用這些組別中的樣本來計算核心菌群"
        )
        if not groups_pick:
            st.warning("請至少選擇一個組別。")
            st.stop()

        # 取子集
        sub_core = merged[merged[label_col].isin(groups_pick)].copy()
        genus_cols = [c for c in sub_core.columns if c not in ["SampleID", label_col]]
        if not genus_cols:
            st.error("找不到菌屬欄位。")
            st.stop()

        # Presence 矩陣（True/False）
        presence_mat = (sub_core[genus_cols] >= presence_thresh)

        # Overall prevalence / mean
        prevalence_overall = presence_mat.mean(axis=0) * 100.0
        mean_overall = sub_core[genus_cols].mean(axis=0, skipna=True)

        # Per-group prevalence（%）
        prev_by_group = (
            presence_mat.join(sub_core[label_col])
            .groupby(label_col)
            .mean() * 100.0
        ).loc[groups_pick]

        # 核心判定
        if mode.startswith("Overall"):
            core_mask = (prevalence_overall >= prev_cut) & (mean_overall >= pooled_mean_min)
        else:
            per_group_ok = (prev_by_group >= prev_cut).all(axis=0)
            core_mask = per_group_ok & (mean_overall >= pooled_mean_min)

        core_list = sorted(prevalence_overall[core_mask].index.tolist())

        # Prevalence–Abundance 散點圖
        st.markdown("### 📈 Prevalence–Abundance（右上角＝核心候選）")
        logy = st.checkbox("Y 軸取 log10（避免高值遮蔽低值）", value=True)

        pa_df = pd.DataFrame({
            "Genus": prevalence_overall.index,
            "Prevalence(%)": prevalence_overall.values,
            "MeanAbundance": mean_overall.values,
            "Core": ["Core" if g in core_list else "Non-core" for g in prevalence_overall.index]
        })
        y_vals = np.log10(pa_df["MeanAbundance"] + 1e-8) if logy else pa_df["MeanAbundance"]
        fig_pa = px.scatter(
            pa_df, x="Prevalence(%)", y=y_vals, color="Core", hover_data=["Genus", "MeanAbundance"],
            color_discrete_map={"Core": "red", "Non-core": "gray"},
            title=f"Presence ≥ {presence_thresh}｜Prevalence ≥ {prev_cut}%｜Pooled mean ≥ {pooled_mean_min}"
        )
        fig_pa.add_vline(x=prev_cut, line_dash="dash")
        if not logy:
            fig_pa.add_hline(y=pooled_mean_min, line_dash="dash")
        st.plotly_chart(fig_pa, use_container_width=True)

        # 各組 Prevalence 熱圖（優先顯示核心菌；若無核心菌則顯示全部）
        st.markdown("### 🔥 各組別 Prevalence 熱圖（%）")
        heat_genus = core_list if len(core_list) > 0 else prevalence_overall.index.tolist()
        heat_df = prev_by_group[heat_genus]
        fig_hm = px.imshow(
            heat_df.values,
            x=heat_genus,
            y=heat_df.index.tolist(),
            aspect="auto",
            labels=dict(x="Genus", y="Group", color="Prevalence (%)"),
            title="Prevalence Heatmap（核心菌優先顯示）"
        )
        st.plotly_chart(fig_hm, use_container_width=True)

        # 核心菌清單 + 下載
        st.markdown(f"### 📜 核心菌清單（共 {len(core_list)} 個）")
        core_table = pd.DataFrame({
            "Genus": core_list,
            "Prevalence_overall(%)": prevalence_overall[core_list].values if len(core_list) > 0 else [],
            "MeanAbundance_overall": mean_overall[core_list].values if len(core_list) > 0 else []
        })
        if not core_table.empty:
            core_table = core_table.sort_values(
                ["Prevalence_overall(%)", "MeanAbundance_overall"], ascending=[False, False]
            )
        st.dataframe(core_table, use_container_width=True)
        st.download_button(
            "📥 下載核心菌清單（CSV）",
            data=core_table.to_csv(index=False).encode(),
            file_name="core_microbiota.csv",
            disabled=core_table.empty
        )

        with st.expander("🔎 偵錯資訊（必要時展開）"):
            st.write("選取組別：", groups_pick)
            st.write("presence_mat 形狀：", presence_mat.shape)
            st.write("prev_by_group（%）前幾列：")
            st.write(prev_by_group.head())

else:
    st.info("請在『📂 上傳資料』分頁上傳檔案。") 
