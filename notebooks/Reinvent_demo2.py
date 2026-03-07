import marimo

__generated_with = "0.18.4"
app = marimo.App()


@app.cell
def _():
    import subprocess
    return (subprocess,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    # _De novo_ Reinforcement Learning with the Reinvent prior

    This is a short demo to
    - Set up a reinforcement  learning run
    - Carry out a reinforcement  learning run
    - Visualize the results with TensorBoard
    - Extract the raw data from TensorBoard and how to work with it
    """)
    return


@app.cell
def _():
    import os
    import shutil

    import pandas as pd
    import seaborn as sns
    import matplotlib.pyplot as plt

    import reinvent
    from reinvent.notebooks import load_tb_data, plot_scalars, get_image, create_mol_grid

    # magic command not supported in marimo; please file an issue to add support
    # %load_ext tensorboard
    return (
        create_mol_grid,
        get_image,
        load_tb_data,
        os,
        pd,
        plot_scalars,
        plt,
        reinvent,
        shutil,
        sns,
    )


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Set up the REINVENT run

    A work directory is defined and created anew (all previous data deleted if the directory already exists).  The prior file is taken from the REINVENT repository and the agent is initially the same as the prior.  A TOML configuration is written out to file in the work directory.
    """)
    return


@app.cell
def _():
    wd = "/tmp/R4_notebooks_output"
    return (wd,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Delete existing working directory (!)

    Skip this step if you have already carried out the RL run but wish to analyse the results further.
    """)
    return


@app.cell
def _(os, shutil, wd):
    shutil.rmtree(wd, ignore_errors=True)
    os.mkdir(wd)
    os.chdir(wd)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Global configuration parameters

    Parameters global to the run:
    - The run type: one of "scoring", "sampling", "transfer_learning", and "staged_learning"
    - The device to run on: "cpu" or "cuda:0" where the number is the device index (needed for ROCm)
    - The output directory for TensorBoard (optional)
    - The configuration file in JSON format (optional)
    """)
    return


@app.cell
def _():
    global_parameters = """
    run_type = "staged_learning"
    device = "mps"   # Use "cpu" if no GPU is available instead "cuda:0"
    tb_logdir = "tb_stage1"
    json_out_config = "_stage1.json"
    """
    return (global_parameters,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Parameters

    Here we specify the model files, the prefix for the output CSV summary file and the batch size for sampling and stochastic gradient descent (SGD).  The batch size is often given in 2^N but there is in now way required.  Typically batch sizes are between 50 and 150.  Batch size effects on SGD and so also the learning rate.  Some experimentation may be required to adjust this but keep in mind that, say, raising the total score as fast as possible is not necessarily the best choice as this may hamper exploration.
    """)
    return


@app.cell
def _(os, reinvent):
    prior_filename = os.path.join(reinvent.__path__[0], "..", "priors", "reinvent.prior")
    agent_filename = prior_filename

    parameters = f"""
    [parameters]

    prior_file = "{prior_filename}"
    agent_file = "{agent_filename}"
    summary_csv_prefix = "stage1"

    batch_size = 100

    use_checkpoint = false
    """
    return (parameters,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Reinforcement Learning strategy
    """)
    return


@app.cell
def _():
    learning_strategy = """
    [learning_strategy]

    type = "dap"
    sigma = 128
    rate = 0.0001
    """
    return (learning_strategy,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ###  Stage setup

    Here we only use a single stage. The aim of this stage is to create an agent which is highly likely to generate "drug-like" molecules (as per QED and Custom Alerts) with no stereocentres

    The stage will terminate when a maximum number of 300 steps is reached.  Termination could occur earlier when the maximum score of 1.0 is exceeded but this is very unlikely to occur.  A checkpoint file is written out which can be used as the agent in a subsequent stage.

    The scoring function is a weighted product of all the scoring components: QED and number of sterecentres.  The latter is used here to avoid stereocentres as they are not supported by the Reinvent prior.  Zero stereocentres aids in downstream 3D task to avoid having to carry out stereocentre enumeration.  Custom alerts is a filter which filters out (scores as zero) all generated compounds which match one of the SMARTS patterns.  Number of sterecentres uses a transformation function to ensure the component score is between 0 and 1.
    """)
    return


@app.cell
def _():
    stages = """
    [[stage]]

    max_score = 1.0
    max_steps = 300

    chkpt_file = 'stage1.chkpt'

    [stage.scoring]
    type = "geometric_mean"

    [[stage.scoring.component]]
    [stage.scoring.component.custom_alerts]

    [[stage.scoring.component.custom_alerts.endpoint]]
    name = "Alerts"

    params.smarts = [
        "[*;r8]",
        "[*;r9]",
        "[*;r10]",
        "[*;r11]",
        "[*;r12]",
        "[*;r13]",
        "[*;r14]",
        "[*;r15]",
        "[*;r16]",
        "[*;r17]",
        "[#8][#8]",
        "[#6;+]",
        "[#16][#16]",
        "[#7;!n][S;!$(S(=O)=O)]",
        "[#7;!n][#7;!n]",
        "C#C",
        "C(=[O,S])[O,S]",
        "[#7;!n][C;!$(C(=[O,N])[N,O])][#16;!s]",
        "[#7;!n][C;!$(C(=[O,N])[N,O])][#7;!n]",
        "[#7;!n][C;!$(C(=[O,N])[N,O])][#8;!o]",
        "[#8;!o][C;!$(C(=[O,N])[N,O])][#16;!s]",
        "[#8;!o][C;!$(C(=[O,N])[N,O])][#8;!o]",
        "[#16;!s][C;!$(C(=[O,N])[N,O])][#16;!s]"
    ]

    [[stage.scoring.component]]
    [stage.scoring.component.QED]

    [[stage.scoring.component.QED.endpoint]]
    name = "QED"
    weight = 0.6


    [[stage.scoring.component]]
    [stage.scoring.component.NumAtomStereoCenters]

    [[stage.scoring.component.NumAtomStereoCenters.endpoint]]
    name = "Stereo"
    weight = 0.4

    transform.type = "left_step"
    transform.low = 0
    """
    return (stages,)


@app.cell
def _(global_parameters, learning_strategy, parameters, stages):
    config = global_parameters + parameters + learning_strategy + stages

    toml_config_filename = "stage1.toml"

    with open(toml_config_filename, "w") as tf:
        tf.write(config)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ##  Start Reinforcement Learning

    This run will take several minutes to finish and timings are provided at the end of the run.  After the run the TensorBoard output can be used to inspect the results of the run
    """)
    return


@app.cell
def _(subprocess):
    # magic command not supported in marimo; please file an issue to add support
    # %%time|
    #! reinvent -l stage1.log $toml_config_filename

    subprocess.call(['reinvent', '-l', 'stage1.log', 'toml_config_filename'])
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Inspect results with TensorBoard

    TensorBoard needs to be started after REINVENT has finished.  Scatter plots of all scoring components are shown (raw and transformed) in the _SCALARS_ tab as well as loss and fractions of valid and duplicate (per batch) SMILES.  The _IMAGES_ tab shows the first sampled molecules for each recorded RL step labelled with the total score for the molecule.

    The loss likelihoods (negative log likelihoods, NLL) for the agent is expected to move away from the prior NLL and also have smaller NLLs than the prior.  This shows that the agent is increasingly producing molecules different from the prior, closer to the signal from the scoring function as requested.  For sample efficiency it is desirable to observed few duplicates and a high number of valid molecules.  Note that results are stochastic and will not be different in every run of RL, even with exactly the same configuration.
    """)
    return


@app.cell
def _():
    # magic command not supported in marimo; please file an issue to add support
    # %tensorboard --bind_all --logdir $wd/tb_stage1_0
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Extract data from TensorBoard

    TensorBoard data can be directly extracted as shown in this section.
    """)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Load the TB data
    """)
    return


@app.cell
def _(load_tb_data, wd):
    ea = load_tb_data(wd)
    return (ea,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Plot all scalars

    All scalar values (except the raw components) are plotted here.  The data is also return as a Pandas dataframe and can so be conveniently used for further analysis or storing to a file.
    """)
    return


@app.cell
def _(ea, plot_scalars):
    df = plot_scalars(ea)
    return (df,)


@app.cell
def _(df):
    df
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Display an image from TB

    Shows (only) the last image from TB.  The image depicts the first 30 molecules generated in the very last step of RL:
    """)
    return


@app.cell
def _(display, ea, get_image):
    img = get_image(ea)
    display(img)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ## Extract data from the CSV file

    The CSV file is generated during the RL run in real time i.e. as soon as the data is available it is written to file. So even if the RL job crashes in the middle of the run some partial data will be available.  The TensorBoard data is a subset of the data in the CSV file.  The CSV file contains in addition all SMILES strings, their state (0=invalid, 1=valid, 2=batch duplicate), the scaffold if a diversity filter has been used.
    """)
    return


@app.cell
def _(os, pd, wd):
    csv_file = os.path.join(wd, 'stage1_1.csv')
    df_1 = pd.read_csv(csv_file)
    df_1
    return (df_1,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Sample efficiency

    Count the number of total invalid and duplicate SMILES and compare to the total number of generated SMILES.
    """)
    return


@app.cell
def _(df_1):
    total_smilies = len(df_1)
    invalids = df_1[df_1['SMILES_state'] == 0]
    total_invalid_smilies = len(invalids)
    duplicates = df_1[df_1['SMILES_state'] == 2]
    total_batch_duplicate_smilies = len(duplicates)
    all_duplicates = df_1[df_1.duplicated(subset=['SMILES'])]
    total_duplicate_smilies = len(all_duplicates)
    print(f'Total number of SMILES generated: {total_smilies}\nTotal number of invalid SMILES generated: {total_invalid_smilies}\nTotal number of batch duplicate SMILES generated: {total_batch_duplicate_smilies}\nTotal number of duplicate SMILES generated: {total_duplicate_smilies}')
    return (all_duplicates,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Display all generated duplicates
    """)
    return


@app.cell
def _(all_duplicates, create_mol_grid, display):
    mol_view = create_mol_grid(all_duplicates)
    display(mol_view)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Display the molecules from the last step
    """)
    return


@app.cell
def _(create_mol_grid, df_1, display):
    last = df_1[df_1['step'] == max(df_1['step'])]
    mol_view_1 = create_mol_grid(last)
    display(mol_view_1)
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""
    ### Plot the NLLs

    The "Target" is the "augmented NLL".
    """)
    return


@app.cell
def _(df_1, plt, sns):
    grouped_df = df_1.groupby('step')
    for label in ('Agent', 'Prior', 'Target'):
        means = grouped_df.aggregate({label: 'mean'})
        X = list(means.index.values)
        sns.scatterplot(means, x=X, y=label, label=label)
    plt.xlabel('step')
    plt.ylabel('NLL')
    plt.show()
    return


@app.cell
def _():
    import marimo as mo
    return (mo,)


if __name__ == "__main__":
    app.run()
