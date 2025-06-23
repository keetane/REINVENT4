import streamlit as st
import os
import datetime # datetimeのインポートを統一
from rdkit import Chem
from rdkit.Chem import Draw, PandasTools
import pandas as pd
import io

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

            # sort options for the selectbox
            desc = ['NLL', 'MW', 'LogP', 'HBD', 'HBA', 'TPSA', 'CtRtBonds', 'CtAmides', 'CtRings', 'CtAromaticRings']
            sort_options=st.sidebar.selectbox(
                "Sort by:",
                options=desc,
                index=0,  # Default to the first option
                key="sort_selector"  # Unique key for the widget
            )
            sort_order = st.sidebar.toggle('Ascending order', key='sort_order', value=True)  # Toggle for ascending/descending order
            # --- メインコンテンツの表示 ---
            # csvの読み込み
            df = pd.read_csv(mols_path, sep=',', encoding='utf-8') if mols_path else pd.DataFrame()
            df = df.sort_values(by=sort_options, ascending=sort_order)

            # filtering parameters
            st.sidebar.markdown("# Filtering Parameters")
            NLL= st.sidebar.select_slider(
                "NLL",
                options=list(range(int(df['NLL'].min()), int(df['NLL'].max()) + 1)),
                value=(int(df['NLL'].min()), int(df['NLL'].max())),
                key="nll_filter"
            )
            MW = st.sidebar.select_slider(
                "MW",
                options=list(range(int(df['MW'].min()), int(df['MW'].max()) + 1)),
                value=(int(df['MW'].min()), int(df['MW'].max())),
                key="mw_filter"
            )
            LogP = st.sidebar.select_slider(
                "LogP",
                options=list(range(int(df['LogP'].min()), int(df['LogP'].max()) + 1)),
                value=(int(df['LogP'].min()), int(df['LogP'].max())),
                key="logp_filter"
            )
            HBD = st.sidebar.select_slider(
                "HBD",
                options=list(range(int(df['HBD'].min()), int(df['HBD'].max()) + 1)),
                value=(int(df['HBD'].min()), int(df['HBD'].max())),
                key="hbd_filter"
            )
            HBA = st.sidebar.select_slider(
                "HBA",
                options=list(range(int(df['HBA'].min()), int(df['HBA'].max()) + 1)),
                value=(int(df['HBA'].min()), int(df['HBA'].max())),
                key="hba_filter"
            )
            TPSA = st.sidebar.select_slider(
                "TPSA",
                options=list(range(int(df['TPSA'].min()), int(df['TPSA'].max()) + 1)),
                value=(int(df['TPSA'].min()), int(df['TPSA'].max())),
                key="tpsa_filter"
            )
            CtRtBonds = st.sidebar.select_slider(
                "CtRtBonds",
                options=list(range(int(df['CtRtBonds'].min()), int(df['CtRtBonds'].max()) + 1)),
                value=(int(df['CtRtBonds'].min()), int(df['CtRtBonds'].max())),
                key="ctrb_filter"
            )
            CtAmides = st.sidebar.select_slider(
                "CtAmides",
                options=list(range(int(df['CtAmides'].min()), int(df['CtAmides'].max()) + 1)),
                value=(int(df['CtAmides'].min()), int(df['CtAmides'].max())),
                key="ctam_filter"
            )
            CtRings = st.sidebar.select_slider(
                "CtRings",
                options=list(range(int(df['CtRings'].min()), int(df['CtRings'].max()) + 1)),
                value=(int(df['CtRings'].min()), int(df['CtRings'].max())),
                key="ctr_filter"
            )
            CtAromaticRings = st.sidebar.select_slider(
                "CtAromaticRings",
                options=list(range(int(df['CtAromaticRings'].min()), int(df['CtAromaticRings'].max()) + 1)),
                value=(int(df['CtAromaticRings'].min()), int(df['CtAromaticRings'].max())),
                key="ctar_filter"
            )
            st.sidebar.text(CtAromaticRings)

            # molfileの生成
            df['ROMol'] = df['SMILES'].apply(lambda x : Chem.MolFromSmiles(x))
            # filtering the dataframe based on the selected parameters
            df = df[(df['NLL'] >= NLL[0]) & (df['NLL'] <= NLL[1])]
            df = df[(df['MW'] >= MW[0]) & (df['MW'] <= MW[1])]
            df = df[(df['LogP'] >= LogP[0]) & (df['LogP'] <= LogP[1])]
            df = df[(df['HBD'] >= HBD[0]) & (df['HBD'] <= HBD[1])]
            df = df[(df['HBA'] >= HBA[0]) & (df['HBA'] <= HBA[1])]
            df = df[(df['TPSA'] >= TPSA[0]) & (df['TPSA'] <= TPSA[1])]
            df = df[(df['CtRtBonds'] >= CtRtBonds[0]) & (df['CtRtBonds'] <= CtRtBonds[1])]
            df = df[(df['CtAmides'] >= CtAmides[0]) & (df['CtAmides'] <= CtAmides[1])]
            df = df[(df['CtRings'] >= CtRings[0]) & (df['CtRings'] <= CtRings[1])]
            df = df[(df['CtAromaticRings'] >= CtAromaticRings[0]) & (df['CtAromaticRings'] <= CtAromaticRings[1])]

            # download button for CSV file
            csv_files = df.drop(columns=['ROMol']).to_csv(index=False, encoding='utf-8')
            st.download_button(
                label="Download CSV",
                data=csv_files,
                file_name=os.path.basename(mols_path) if mols_path else "molecules.csv",
                mime="text/csv",
                key="download_csv"
            )
            # st.dataframe(df) if mols_path else st.info("CSVファイルが選択されていません。")

            img = Draw.MolsToGridImage(
                df['ROMol'],
                molsPerRow=3,
                subImgSize=(400, 300),
                legends=df['NLL'].astype(str).tolist() if 'NLL' in df.columns else None
            )

            st.image(img, caption="Molecules from the selected CSV file")

        else:
            st.sidebar.info("CSVファイルをサイドバーから選択してください。")
    else:
        st.sidebar.warning("resultsディレクトリ内にCSVファイルが見つかりません。")
else:
    st.sidebar.error("resultsディレクトリが見つからないか、アクセスできません。")
    st.error("resultsディレクトリが見つからないか、アクセスできません。")


