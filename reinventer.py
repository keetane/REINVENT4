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
            desc = ['NLL', 'MW', 'LogP', 'HBD', 'HBA', 'TPSA', 'CtRotB', 'CtAmides', 'CtRings', 'CtAromaticRings', 'Tanimoto_Similarity_to_Parent']
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

            filter_values = {}

            # Define which columns are float and which are int
            float_columns = ['NLL', 'MW', 'LogP', 'TPSA']
            int_columns = ['HBD', 'HBA', 'CtRotBonds', 'CtAmides', 'CtRings', 'CtAromaticRings']

            # Handle float columns with st.slider
            for col in float_columns:
                if col in df.columns and not df[col].isnull().all():
                    min_val, max_val = float(df[col].min()), float(df[col].max())
                    if min_val >= max_val:
                        st.sidebar.text(f"{col}: {min_val:.2f} (single value)")
                        filter_values[col] = (min_val, max_val)
                    else:
                        filter_values[col] = st.sidebar.slider(
                            col,
                            min_value=min_val,
                            max_value=max_val,
                            value=(min_val, max_val),
                            key=f"{col}_filter"
                        )

            # Handle int columns with st.select_slider
            for col in int_columns:
                if col in df.columns and not df[col].isnull().all():
                    min_val, max_val = int(df[col].min()), int(df[col].max())
                    if min_val >= max_val:
                        st.sidebar.text(f"{col}: {min_val} (single value)")
                        filter_values[col] = (min_val, max_val)
                    else:
                        options = list(range(min_val, max_val + 1))
                        filter_values[col] = st.sidebar.select_slider(
                            col,
                            options=options,
                            value=(min_val, max_val),
                            key=f"{col}_filter"
                        )
            
            # Special handling for Tanimoto
            if 'Tanimoto_Similarity_to_Parent' in df.columns and not df['Tanimoto_Similarity_to_Parent'].isnull().all():
                col = 'Tanimoto_Similarity_to_Parent'
                min_val, max_val = float(df[col].min()), float(df[col].max())
                if min_val >= max_val:
                    st.sidebar.text(f"Tanimoto Similarity: {min_val:.2f} (single value)")
                    filter_values[col] = (min_val, max_val)
                else:
                    filter_values[col] = st.sidebar.slider(
                        "Tanimoto Similarity to Parent",
                        min_value=min_val,
                        max_value=max_val,
                        value=(min_val, max_val),
                        step=0.01,
                        key="tanimoto_filter"
                    )

            # molfileの生成
            df['ROMol'] = df['SMILES'].apply(lambda x : Chem.MolFromSmiles(x))
            
            # filtering the dataframe based on the selected parameters
            for col, val_range in filter_values.items():
                if val_range:
                    df = df[(df[col] >= val_range[0]) & (df[col] <= val_range[1])]

            # --- メインコンテンツの表示 ---
            col1, col2, col3 = st.columns(3)

            with col1:
                st.metric(label="Molecules", value=len(df))

            with col2:
                csv_files = df.drop(columns=['ROMol']).to_csv(index=False, encoding='utf-8')
                st.download_button(
                    label="Download CSV",
                    data=csv_files,
                    file_name=os.path.basename(mols_path) if mols_path else "molecules.csv",
                    mime="text/csv",
                    key="download_csv"
                )

            with col3:
                # SDFをストリームに書き出す
                sdf_buffer = io.StringIO()
                PandasTools.WriteSDF(df, sdf_buffer, molColName='ROMol', properties=df.columns.drop('ROMol').tolist())
                st.download_button(
                    label="Download SDF",
                    data=sdf_buffer.getvalue().encode("utf-8"),
                    file_name=os.path.basename(mols_path).replace('.csv', '.sdf') if mols_path else "molecules.sdf",
                    mime="chemical/x-mdl-sdfile",
                    key="download_sdf"
                )
            # st.dataframe(df) if mols_path else st.info("CSVファイルが選択されていません。")

            img = Draw.MolsToGridImage(
                df['ROMol'].head(100),
                molsPerRow=3,
                subImgSize=(400, 300),
                legends=df[sort_options].head(100).astype(str).tolist() if sort_options in df.columns else None
            )

            st.image(img, caption="Molecules from the selected CSV file")

        else:
            st.sidebar.info("CSVファイルをサイドバーから選択してください。")
    else:
        st.sidebar.warning("resultsディレクトリ内にCSVファイルが見つかりません。")
else:
    st.sidebar.error("resultsディレクトリが見つからないか、アクセスできません。")
    st.error("resultsディレクトリが見つからないか、アクセスできません。")


