
# Volcano Analyzer V10

Two-group Volcano plot app with FDR/p-value toggle, Top5, boxplot, and core microbiota analysis.

## 🚀 Deploy on Streamlit Community Cloud
1. Push this repo to GitHub.
2. Go to https://share.streamlit.io/ → **New app** → Select this repo.
3. Main file: `app.py` → Deploy.

## 🖥️ Run locally
```bash
python -m venv .venv && . .venv/bin/activate  # Windows: .venv\Scripts\activate
pip install -r requirements.txt
streamlit run app.py
```

## 📂 Data format
- **Genus abundance (.tsv)**: first column `Genus`, following columns are sample IDs.
- **Metadata (.xlsx)**: include a sample-ID column (select it in the UI) and a grouping column.

MIT © 2025
