import os
import subprocess
import pandas as pd
from rdkit.Chem import PandasTools, Draw
from rdkit import Chem
import streamlit as st
from datetime import datetime
from rdkit.Chem import Recap
import pubchempy as pcp
import datetime

# Get the current timestamp
now = datetime.datetime.now()
# Format the timestamp as a string
time = now.strftime("%Y%m%d_%H%M")

st.set_page_config(
    page_title="REINVENTer 4 Drug Discovery",
    page_icon=":pill:",
    layout="wide",
    initial_sidebar_state="expanded",
)


# Set up directories relative to the script's location
wd = st.session_state.wd  # Use the working directory from session state
os.makedirs(wd, exist_ok=True)  # Create directory if it doesn't exist
reinvent_dir = os.getcwd()  # Adjusted REINVENT4 directory path
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

# Set up file paths for priors
priors_dir = os.path.join(reinvent_dir, "priors")
reinvent = os.path.join(priors_dir, "reinvent.prior")
lib = os.path.join(priors_dir, "libinvent.prior")
link = os.path.join(priors_dir, "linkinvent.prior")
mol2mol_high = os.path.join(priors_dir, "mol2mol_high_similarity.prior")
mol2mol_med = os.path.join(priors_dir, "mol2mol_medium_similarity.prior")
mol2mol_mmp = os.path.join(priors_dir, "mol2mol_mmp.prior")
mol2mol_scaffold_generic = os.path.join(priors_dir, "mol2mol_scaffold_generic.prior")
mol2mol_scaffold = os.path.join(priors_dir, "mol2mol_scaffold.prior")
mol2mol_similarity = os.path.join(priors_dir, "mol2mol_similarity.prior")
pubchem = os.path.join(priors_dir, "pubchem_ecfp4_with_count_with_rank_reinvent4_dict_voc.prior")

# Fetch SMILES from PubChem
st.header("Parent Molecule Drawer")

# separate with columns
col1, col2 = st.columns([2, 1])
with col1:

    compound_name = st.text_input("Enter compound name", value="ruxolitinib")
    if st.button("Fetch SMILES from PubChem into text_area"):
        try:
            smiles = pcp.get_compounds(compound_name, 'name')[0].isomeric_smiles
        except Exception as e:
            st.error(f"Error fetching SMILES for {compound_name}: {e}")
    else:
        smiles = "C1CCC(C1)[C@@H](CC#N)N2C=C(C=N2)C3=C4C=CNC4=NC=N3"  # Default blank

    # Display fetched SMILES
    smiles = st.text_area("enter SMILES", value=smiles, height=68)

with col2:
    st.image(Draw.MolToImage(Chem.MolFromSmiles(smiles), size=(300, 300)), caption="Parent Molecule")

# save SMILES of parent molecule from text_area
with open(f"{input_dir}/parent.smi", "w") as f:
    f.write(f"{smiles}\n")

# Recap decomposition and molecule visualization
st.header("Children by Recap Decomposition")

mol = Chem.MolFromSmiles(smiles)
recap_tree = Recap.RecapDecompose(mol)
warhead_mols = []
warhead_smiless = []
child_mols = []
child_smiless = []

# ダミー原子にatom map番号を付ける関数
def relabel_dummy_atoms(mol, map_num=1):
    """mol内のすべてのダミー原子(*)にatom map番号を付ける → [*:1]形式に"""
    rw_mol = Chem.RWMol(mol)
    for atom in rw_mol.GetAtoms():
        if atom.GetAtomicNum() == 0:  # dummy atom (*)
            atom.SetAtomMapNum(map_num)
    return rw_mol.GetMol()

# 分子を1stepだけ分解し、ラベルを変換
for smiles, node in recap_tree.children.items():
    mol = node.mol
    warhead_mols.append(mol)
    warhead_smiless.append(Chem.MolToSmiles(mol))
    modified = relabel_dummy_atoms(mol, map_num=1)
    child_mols.append(modified)
    child_smiless.append(Chem.MolToSmiles(modified))

BB_namelist = [f'BB{i+1}' for i in range(len(warhead_mols))]
st.image(Draw.MolsToGridImage(warhead_mols, molsPerRow=3, subImgSize=(600,300), legends=BB_namelist))


# Sidebar for input parameters
st.sidebar.header("Sampling Parameters")
num_mols = st.sidebar.number_input("Number of SMILES", min_value=1, value=155)
device = st.sidebar.selectbox("Device", ["cpu", "cuda", "mps"])  # "cpu"を先頭に

unique_molecules = st.sidebar.checkbox("Unique Molecules", value=True)
randomize_smiles = st.sidebar.checkbox("Randomize SMILES", value=True)
overwrite = st.sidebar.checkbox("Overwrite", value=True)

# Allow users to select a model file from priors_dir
model_files = {
    "High Similarity": mol2mol_high,
    "Medium Similarity": mol2mol_med,
    "MMP": mol2mol_mmp,
    "Scaffold Generic": mol2mol_scaffold_generic,
    "Scaffold": mol2mol_scaffold,
    "Similarity": mol2mol_similarity,
    "PubChem": pubchem,
    "Reinvent": reinvent,
}
model_file = [] #st.sidebar.selectbox("Model File", options=model_files.keys(), format_func=lambda x: x)
selected_model_file = model_files['Reinvent']

st.sidebar.header('mol2mol options')
# other options for mol2mol sampling
sample_stategies = st.sidebar.selectbox("Sampling Strategy", ["beamsearch", 'multinomial'])
temperatures = st.sidebar.number_input("Temperature", min_value=0.0, max_value=1.0, value=1.0, step=0.1)

# # Generate TOML file
# def generate_toml(method=selected_model_file, num_smiles=num_mols, smiles_file=None, device=device, show=False, sample_strategy=None, temperature=1.0):
#     if method == link:
#         filename = 'LinkInvent'
#     elif method == lib:
#         filename = 'LibInvent'
#     else:
#         filename = method.split("/")[-1].replace(" ", "_")[:-6]
        
#     if overwrite is not True:
#         filename = method.split("/")[-1].replace(" ", "_") + '_' + time
#     else:
#         pass
#     output_file = os.path.join(results_dir, f"{filename}.csv")
#     toml_content = f"""
# run_type = "sampling"
# device = "{device}"
# json_out_config = "{sampling_log}/_sampling.json"

# [parameters]
# model_file = "{method}"
# output_file = "{output_file}"
# num_smiles = {num_smiles}
# """
#     if unique_molecules is True:
#         toml_content += f'''
#         unique_molecules = true
#         '''
#     if randomize_smiles is True:
#         toml_content += f'''
#         randomize_smiles = true
#         '''
#     if smiles_file is not None:
#         toml_content += f'''
#         smiles_file = "{smiles_file}"  # 1 compound per line
#         '''
#     # Fix: use toml_content instead of undefined toml
#     if sample_strategy == 'beamsearch':
#         toml_content += f'''
#         sample_strategy = "beamsearch"  # multinomial or beamsearch (deterministic)
#         '''
#     elif sample_strategy == 'multinomial':
#         toml_content += f'''
#         sample_strategy = "multinomial"  # multinomial or beamsearch (deterministic)
#         temperature = {temperature} # temperature in multinomial sampling
#         '''

#     toml_path = os.path.join(toml_dir, "sampling.toml")
#     with open(toml_path, "w") as f:
#         f.write(toml_content)
#     return toml_path, output_file

# # Run REINVENT4
# def run_reinvent(toml_path):
#     log_file = os.path.join(sampling_log, "sampling.log")
#     subprocess.call(["reinvent", "-l", log_file, toml_path])

# Generate TOML and run REINVENT4
def reinvent(method=selected_model_file, num_smiles=num_mols, smiles_file=None, device=device, show=False, sample_strategy=None, temperature=1.0):
    if method == link:
        filename = 'LinkInvent'
    elif method == lib:
        filename = 'LibInvent'
    else:
        filename = method.split("/")[-1].replace(" ", "_")[:-6]
        
    if overwrite is not True:
        filename = method.split("/")[-1].replace(" ", "_") + '_' + time
    else:
        pass
    output_file = os.path.join(results_dir, f"{filename}.csv")
    toml_content = f"""
run_type = "sampling"
device = "{device}"
json_out_config = "{sampling_log}/_sampling.json"

[parameters]
model_file = "{method}"
output_file = "{output_file}"
num_smiles = {num_smiles}
"""
    if unique_molecules is True:
        toml_content += f'''
        unique_molecules = true
        '''
    if randomize_smiles is True:
        toml_content += f'''
        randomize_smiles = true
        '''
    if smiles_file is not None:
        toml_content += f'''
        smiles_file = "{smiles_file}"  # 1 compound per line
        '''
    # Fix: use toml_content instead of undefined toml
    if sample_strategy == 'beamsearch':
        toml_content += f'''
        sample_strategy = "beamsearch"  # multinomial or beamsearch (deterministic)
        '''
    elif sample_strategy == 'multinomial':
        toml_content += f'''
        sample_strategy = "multinomial"  # multinomial or beamsearch (deterministic)
        temperature = {temperature} # temperature in multinomial sampling
        '''

    toml_path = os.path.join(toml_dir, "sampling.toml")
    with open(toml_path, "w") as f:
        f.write(toml_content)
    # return toml_path, output_file

    # Run REINVENT4
    log_file = os.path.join(sampling_log, "sampling.log")
    subprocess.call(["reinvent", "-l", log_file, toml_path])
    df = pd.read_csv(output_file, sep=',')
    df['ROMol']= df['SMILES'].map(lambda x: Chem.MolFromSmiles(x))
    # 分子のdescripterを追加
    df['MW'] = df['ROMol'].map(lambda x: Chem.Descriptors.MolWt(x) if x else None)
    df['LogP'] = df['ROMol'].map(lambda x: Chem.Descriptors.MolLogP(x) if x else None)
    df['HBD'] = df['ROMol'].map(lambda x: Chem.Descriptors.NumHDonors(x) if x else None)
    df['HBA'] = df['ROMol'].map(lambda x: Chem.Descriptors.NumHAcceptors(x) if x else None)
    df['TPSA'] = df['ROMol'].map(lambda x: Chem.Descriptors.TPSA(x) if x else None)
    df['CtRtbonds'] = df['ROMol'].map(lambda x: Chem.rdMolDescriptors.CalcNumRotatableBonds(x) if x else None)
    df['CtRings'] = df['ROMol'].map(lambda x: Chem.rdMolDescriptors.CalcNumRings(x) if x else None)
    df['CtAmides'] = df['ROMol'].map(lambda x: Chem.rdMolDescriptors.CalcNumAmideBonds(x) if x else None)
    df['CtAromaticRings'] = df['ROMol'].map(lambda x: Chem.rdMolDescriptors.CalcNumAromaticRings(x) if x else None)
    df['CtHeteroatoms'] = df['ROMol'].map(lambda x: Chem.rdMolDescriptors.CalcNumHeteroatoms(x) if x else None)
    df.drop(columns=['ROMol'], inplace=True)  # Remove the ROMol column if not needed
    df.to_csv(output_file, index=False)  # Save the updated DataFrame to CSV

# 3columns for sampling
st.header("ReInvent Sampling")
col1, col2, col3 = st.columns([1, 1, 1])

# mol2mol sampling
with col1:
    st.markdown('#### Mol2Mol Sampling')
    model_file = st.selectbox("Model File", options=model_files.keys(), format_func=lambda x: x)
    selected_model_file = model_files[model_file]
    if st.button('Mol2Mol Sampling from Parent', key="mol2mol_sampling_parent_btn"):
        # Generate TOML file for LibInvent
        reinvent(
            method=selected_model_file,
            smiles_file=f'{input_dir}/parent.smi',
            num_smiles=num_mols,
            device=device,
            sample_strategy=sample_stategies,
            temperature=temperatures
        )
        st.success("Sampling completed!")

# LibInvent用の親分子を選択
with col2:
    st.markdown('#### LibInvent Sampling')
    selected_child = st.multiselect(
        "Select Child Molecule", options=BB_namelist, default=["BB2"], key="child_multiselect"
    )
    if selected_child:
        selected_child_index = BB_namelist.index(selected_child[0])
        selected_child_smiles = child_smiless[selected_child_index]
        with open(f"{input_dir}/child.smi", "w") as f:
            f.write(f"{selected_child_smiles}\n")
    else:
        selected_child_smiles = ""
        st.text("No parent selected.")

    if st.button('LibInvent Sampling', key="libinvent_sampling_btn"):
        # Generate TOML file for LibInvent
        reinvent(
            method=lib,
            smiles_file=f'{input_dir}/child.smi',
            num_smiles=num_mols,
            device=device
        )
        st.success("Sampling completed!")

# LinkInvent用のwarheadを選択
with col3:
    st.markdown('#### LinkInvent')
    selected_warhead = st.multiselect(
        "Select Warhead", options=BB_namelist, default=["BB1", "BB2"], key="warhead_multiselect"
    )
    if selected_warhead:
        selected_warhead_indices = [BB_namelist.index(w) for w in selected_warhead]
        selected_warhead_smiles_list = [warhead_smiless[i] for i in selected_warhead_indices]
        with open(f"{input_dir}/warheads.smi", "w") as f:
            f.write('|'.join(selected_warhead_smiles_list))
    else:
        selected_warhead_smiles = ""
        st.text("No warhead selected.")

    if st.button('LinkInvent Sampling', key="linkinvent_sampling_btn"):
        # Generate TOML file for LibInvent
        reinvent(
            method=link,
            smiles_file=f'{input_dir}/warheads.smi',
            num_smiles=num_mols,
            device=device
        )
        st.success("Sampling completed!")