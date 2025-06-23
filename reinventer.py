import streamlit as st
import os
import datetime # datetimeのインポートを統一
from rdkit import Chem
from rdkit.Chem import Draw, PandasTools
import pandas as pd

st.set_page_config(
    page_title="REINVENTer 4 Drug Discovery",
    page_icon=":pill:",
    layout="wide",
    initial_sidebar_state="expanded",
)

# --- サイドバーの設定 ---
st.sidebar.header("Re:Inventer 4 Drug Discovery")
st.sidebar.markdown("[Learn more about REINVENT4](https://jcheminf.biomedcentral.com/articles/10.1186/s13321-024-00812-5)")
st.sidebar.title("Working Directory")

# current directoryの取得
current_dir = os.getcwd()

# current directoryの下に1_wdディレクトリが存在するか確認
if not os.path.exists(os.path.join(current_dir, "1_wd")):
    os.makedirs(os.path.join(current_dir, "1_wd"))
st.session_state.wd = os.path.join(current_dir, "1_wd")

# セッションステートに 'user_subfolder' がなければ初期値を 'all_users' に設定
if 'user_subfolder' not in st.session_state:
    st.session_state.user_subfolder = "all_users"

# サイドバーにサブフォルダ名を入力するテキストボックス
user_subfolder_input = st.sidebar.text_input(
    st.session_state.wd,
    value=st.session_state.user_subfolder
)

# 入力値が変更されたらセッションステートを更新し、再実行
if user_subfolder_input != st.session_state.user_subfolder:
    st.session_state.user_subfolder = user_subfolder_input
    st.rerun()

# 1_wdディレクトリの下にユーザーのサブフォルダを作成
if 'all_users' not in st.session_state:
    # ユーザーのサブフォルダ名が空でなければ、それを使用
    if st.session_state.user_subfolder:
        st.session_state.wd = os.path.join(current_dir, "1_wd", st.session_state.user_subfolder)
    else:
        st.session_state.wd = os.path.join(current_dir, "1_wd", "all_users")

# # --- ディレクトリパスの構築と作成 ---
# # ユーザーのホームディレクトリを展開し、REINVENT4/wdまでのベースパスを構築
# base_path = os.path.expanduser(os.path.join("~", "Documents", "apps", "REINVENT4", "wd"))

# # 最終的な作業ディレクトリパスを決定
# # もしユーザーがサブフォルダ名を入力していればそれを結合、そうでなければベースパスのみ
# st.session_state.wd = os.path.join(base_path, st.session_state.user_subfolder)

# ディレクトリを作成 (既に存在すれば何もしない)
wd = st.session_state.wd  # Use the working directory from session state
os.makedirs(wd, exist_ok=True)
script_dir = os.path.dirname(os.path.abspath(__file__))  # Get the script's directory
input_dir = os.path.join(wd, 'input')  # Input directory path
os.makedirs(input_dir, exist_ok=True)  # Create input directory if it doesn't exist
model_dir = os.path.join(wd, 'model')  # Model directory path
os.makedirs(model_dir, exist_ok=True)  # Create model directory if it doesn't exist
results_dir = os.path.join(wd, 'results')  # Results directory path
os.makedirs(results_dir, exist_ok=True)  # Create results directory if it doesn't exist
toml_dir = os.path.join(wd, 'toml')  # TOML directory path
os.makedirs(toml_dir, exist_ok=True)  # Create TOML directory if it doesn't exist
sampling_log = os.path.join(results_dir, 'log')  # Sampling log directory path
os.makedirs(sampling_log, exist_ok=True)  # Create sampling log directory if it doesn't exist

# --- CSVファイル選択のためのselectbox ---
st.sidebar.markdown("### Results CSV File Selection")

results_dir_to_search = os.path.join(st.session_state.wd, 'results')
csv_files = []
mols_path = None # Initialize mols_path

# Debugging: Show the directory being searched

if os.path.exists(results_dir_to_search) and os.path.isdir(results_dir_to_search):
    for root, _, files in os.walk(results_dir_to_search):
        for f in files:
            if f.endswith('.csv'):
                csv_files.append(os.path.join(root, f))
    
    csv_files.sort() # Sort files for better user experience

    if csv_files:
        # Create a list of relative paths for display in the selectbox
        display_files = [os.path.relpath(f, results_dir_to_search) for f in csv_files]
        display_files.insert(0, "--- ファイルを選択してください ---") # Add a default "Select a file" option
        
        selected_file = st.sidebar.selectbox(
            "Choose a CSV file:", 
            options=display_files,
            key="csv_selector" # Unique key for the widget
        )

        if selected_file != "--- ファイルを選択してください ---":
            # Get the absolute path of the selected file
            mols_path_index = display_files.index(selected_file) - 1 # Adjust index because of the inserted default option
            mols_path = csv_files[mols_path_index]
            # Debugging: Show the assigned mols_path

            # --- メインコンテンツの表示 ---
            df = pd.read_csv(mols_path, sep=',', encoding='utf-8') if mols_path else pd.DataFrame()
            df = df.sort_values(by='NLL', ascending=True)
            df['ROMol'] = df['SMILES'].apply(lambda x : Chem.MolFromSmiles(x))
            # st.dataframe(df) if mols_path else st.info("CSVファイルが選択されていません。")
            img = Draw.MolsToGridImage(
                df['ROMol'],
                molsPerRow=3,
                subImgSize=(300, 200),
                legends=df['NLL'].astype(str).tolist() if 'NLL' in df.columns else None
            )

            st.image(img, caption="Molecules from the selected CSV file", use_column_width=True)

        else:
            st.sidebar.info("CSVファイルをサイドバーから選択してください。")
    else:
        st.sidebar.warning("resultsディレクトリ内にCSVファイルが見つかりません。")
else:
    st.sidebar.error("resultsディレクトリが見つからないか、アクセスできません。")
    st.error("resultsディレクトリが見つからないか、アクセスできません。")



