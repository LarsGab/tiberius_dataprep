import os
import glob
import re
import scipy

from matplotlib.gridspec import GridSpec
import matplotlib.pyplot as plt
def read_training_stats(training_dir, species_list=[]):
    """
    Reads stats for all epochs for a given training number.
    Assumes folder structure: train<number>/epoch_<numb>.stats
    File format example:
    
    #-----------------| Sensitivity | Precision  |
            Base level:    88.0     |    87.3    |
            Exon level:    68.8     |    64.9    |
          Intron level:    74.2     |    70.3    |
    Intron chain level:    33.0     |    32.8    |
      Transcript level:    38.6     |    36.4    |
           Locus level:    41.7     |    36.4    |
    
    Returns:
        dict: Mapping epoch numbers (int) to a dict mapping each level name to its metrics.
              For example:
              {
                  0: {
                      "Base": {"sensitivity": 88.0, "precision": 87.3},
                      "Exon": {"sensitivity": 68.8, "precision": 64.9},
                      ...
                  },
                  1: { ... },
                  ...
              }
    """
    directory = training_dir
    
    epoch_stats = {}
    for s in species_list:
        epoch_stats[s] = {}
        stats_files = sorted(glob.glob(os.path.join(directory, f"*.epoch_*.{s}.stats")))  
        # Regular expression to capture the level name, sensitivity, and precision.
        # It looks for lines like:
        # "        Base level:    88.0     |    87.3    |"
        pattern = re.compile(r'^\s*(.*?)\s+level:\s+([\d.]+)\s+\|\s+([\d.]+)', re.IGNORECASE)
        for file_path in stats_files:
            # Extract epoch number from the filename: "epoch_<number>.stats"
            match_epoch = re.search(r'(\d+).epoch_(\d+)\.', file_path)
            if not match_epoch:
                continue
            train_num = int(match_epoch.group(1))
            epoch_num = int(match_epoch.group(2))
            
            level_metrics = {}
            with open(file_path, 'r') as f:
                for line in f:
                    m = pattern.match(line)
                    if m:
                        level_name = m.group(1).strip()
                        sensitivity = float(m.group(2))
                        precision = float(m.group(3))
                        level_metrics[level_name] = {"sensitivity": sensitivity, "precision": precision}
            if train_num not in epoch_stats[s]:
                epoch_stats[s][train_num] = {}            
            epoch_stats[s][train_num][epoch_num] = level_metrics

    return epoch_stats

def get_best_epoch(stats, metrics=["Transcript", "Exon"], strict=True):
    curr_best = -1
    epoch = [-1,-1]
    for i in range(500):
        if i not in stats[0].keys():
            break 
        for j in range(500):
            curr = []
            for s in stats:
                if j not in s[i].keys():
                    curr.append(-500)
                else:
                    curr.append(sum([s[i][j][m]["sensitivity"] + s[i][j][m]["precision"] for m in metrics]))
            av = sum(curr)/len(stats)
            epoch = [i,j] if av > curr_best else epoch
            curr_best = av if av > curr_best else curr_best
    return epoch, curr_best
            

def plot_training_metrics(stats, title="", out=""):
    """
    Plots sensitivity and precision for each level over epochs.

    Parameters:
        stats (dict): A dictionary where each key is an epoch (int) and each value is a dictionary
                      that maps level names (str) to their metrics (a dict with keys 'sensitivity' and 'precision').

    The function creates a figure with two subplots:
      - The top subplot plots sensitivity values over epochs.
      - The bottom subplot plots precision values over epochs.
    Each level is plotted with its own line for easy comparison.
    """
    # Get sorted list of epochs
    
    epoch_stats_sorted = {}
    for i in range(len(stats.keys())):
        for ep in sorted(stats[i].keys()):
            epoch_stats_sorted[len(epoch_stats_sorted.keys())] = stats[i][ep]
    stats = epoch_stats_sorted

    epochs = sorted(stats.keys())
        
    # Create figure with two subplots
    fig, axs = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    if title:
        fig.suptitle(title, fontstyle="italic")
    # Plot sensitivity for each level
    for level in ["Base", "Exon", "Intron", "Transcript"]:
        sens_values = []
        for ep in epochs:
            # If the level isn't found for a given epoch, use None (or you could use np.nan)
            metrics = stats[ep].get(level, {})
            sens_values.append(metrics.get("sensitivity", None))
        axs[0].plot(epochs, sens_values,  label=level)
    
    axs[0].set_ylabel("Sensitivity [%]")
    axs[0].set_title("Sensitivity")
    axs[0].grid(True)
    axs[0].legend(loc="best")
    
    # Plot precision for each level
    for level in ["Base", "Exon", "Intron", "Transcript"]:
        prec_values = []
        for ep in epochs:
            metrics = stats[ep].get(level, {})
            prec_values.append(metrics.get("precision", None))
        axs[1].plot(epochs, prec_values,  label=level)
    
    axs[1].set_xlabel("Epoch")
    axs[1].set_ylabel("Precision [%]")
    axs[1].set_title("Precision")
    axs[1].grid(True)
    axs[1].legend(loc="best")
    
    for i in range(2):
        axs[i].set_ylim((0,100))
    plt.tight_layout()
    if out:
        plt.savefig(out, dpi=300)
    plt.show()


def plot_accuracies(stats_list, run_labels=None):
    """
    Plots the accuracies (computed as (sensitivity + precision)/2) on Transcript and Exon levels,
    for a list of stat dictionaries from different training runs.

    Each stats dictionary should be structured as:
      {
         epoch_number: {
             "Transcript": {"sensitivity": value, "precision": value},
             "Exon": {"sensitivity": value, "precision": value},
             ...
         },
         ...
      }
    
    Parameters:
      stats_list (list): A list of stat dictionaries, each corresponding to a training run.
      run_labels (list, optional): A list of labels for the training runs. If not provided,
                                   labels "Run 1", "Run 2", ... will be used.

    The function creates a figure with two subplots:
      - The top subplot shows Transcript level accuracies over epochs.
      - The bottom subplot shows Exon level accuracies over epochs.
    """
    # Create default labels if none are provided
    if run_labels is None:
        run_labels = [f"Run {i+1}" for i in range(len(stats_list))]
    
    # Create two subplots sharing the same x-axis
    fig, axs = plt.subplots(2, 1, figsize=(10, 8), sharex=True)
    
    # Iterate over each stats dictionary (i.e. each training run)
    for stats, label in zip(stats_list, run_labels):
        # Get sorted list of epochs
        epochs = sorted(stats.keys())
        
        transcript_accuracies = []
        exon_accuracies = []
        
        for ep in epochs:
            metrics = stats[ep]
            # Compute Transcript accuracy if available, else use None
            if "Transcript" in metrics:
                trans = metrics["Transcript"]
                transcript_accuracies.append(scipy.stats.hmean([trans["sensitivity"], trans["precision"]]))
            else:
                transcript_accuracies.append(None)
            
            # Compute Exon accuracy if available, else use None
            if "Exon" in metrics:
                exon = metrics["Exon"]
                exon_accuracies.append(scipy.stats.hmean([exon["sensitivity"] + exon["precision"]]))
            else:
                exon_accuracies.append(None)
        
        axs[0].plot(epochs[:20], transcript_accuracies[:20],  label=label)
        axs[1].plot(epochs[:20], exon_accuracies[:20],  label=label)
    
    # Configure the Transcript subplot
    axs[0].set_ylabel("Transcript Accuracy [%]")
    axs[0].set_title("Transcript Level Accuracy over Epochs")
    axs[0].grid(True)
    axs[0].legend(loc="best")
    
    # Configure the Exon subplot
    axs[1].set_xlabel("Epoch")
    axs[1].set_ylabel("Exon Accuracy [%]")
    axs[1].set_title("Exon Level Accuracy over Epochs")
    axs[1].grid(True)
    axs[1].legend(loc="best")
    
    plt.tight_layout()
    plt.show()

def plot_gffcompare_accuracy(file_paths, labels=None):
    """
    Given a list of file paths to gffcompare summary outputs, extract the
    Exon-level and Transcript-level Sensitivity and Precision and plot them.

    Parameters
    ----------
    file_paths : list of str
        Paths to files containing gffcompare summaries.
    labels : list of str, optional
        Labels for each file (must be same length as file_paths). If None,
        the basename of each file will be used.

    Returns
    -------
    None
        Displays a Matplotlib figure with two subplots.
    """
    # prepare labels
    if labels is None:
        labels = [os.path.basename(p) for p in file_paths]
    elif len(labels) != len(file_paths):
        raise ValueError("`labels` must be same length as `file_paths`")

    # storage for metrics
    exon_sens, exon_prec = [], []
    tx_sens, tx_prec = [], []

    # regex patterns
    exon_pat = re.compile(r'Exon level:\s*([\d\.]+)\s*\|\s*([\d\.]+)')
    tx_pat   = re.compile(r'Transcript level:\s*([\d\.]+)\s*\|\s*([\d\.]+)')

    # parse each file
    for path in file_paths:
        with open(path) as f:
            text = f.read()

        m_exon = exon_pat.search(text)
        m_tx   = tx_pat.search(text)

        if not m_exon or not m_tx:
            raise RuntimeError(f"Could not find Exon/Transcript lines in {path!r}")

        # convert to floats
        sens_exon = float(m_exon.group(1))
        prec_exon = float(m_exon.group(2))
        sens_tx   = float(m_tx.group(1))
        prec_tx   = float(m_tx.group(2))

        exon_sens.append(sens_exon)
        exon_prec.append(prec_exon)
        tx_sens.append(sens_tx)
        tx_prec.append(prec_tx)

    # plotting
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 6))

    # Exon plot
    ax1.scatter(exon_prec, exon_sens)
    for x, y, lbl in zip(exon_prec, exon_sens, labels):
        ax1.annotate(lbl, (x, y), textcoords="offset points", xytext=(-5,-5))
    ax1.set_xlabel("Precision [%]")
    ax1.set_ylabel("Sensitivity [%]")
    ax1.set_title("Exon Accuracy")
    ax1.set_ylim((0,100))
    ax1.set_xlim((0,100))

    # Transcript plot
    ax2.scatter(tx_prec, tx_sens)
    for x, y, lbl in zip(tx_prec, tx_sens, labels):
        ax2.annotate(lbl, (x, y), textcoords="offset points", xytext=(-5,-5))
    ax2.set_xlabel("Precision [%]")
    ax2.set_ylabel("Sensitivity [%]")
    ax2.set_title("Transcript Accuracy")
    ax2.set_ylim((0,100))
    ax2.set_xlim((0,100))

    plt.tight_layout()
    plt.show()

def draw_metrics_table(ax, labels, exon_sens, exon_prec, exon_f1, tx_sens, tx_prec, tx_f1):
    """
    Draws a table of Sensitivity, Precision, and F1 metrics on the given Axes.
    """
    ax.axis("off")
    # prepare table content
    cell_text = [
        [labels[i],
         f"{exon_sens[i]:.1f}",
         f"{exon_prec[i]:.1f}",
         f"{exon_f1[i]:.1f}",
         f"{tx_sens[i]:.1f}",
         f"{tx_prec[i]:.1f}",
         f"{tx_f1[i]:.1f}"]
        for i in range(len(labels))
    ]
    col_labels = [
        "",
        "Exon Sens [%]", "Exon Prec [%]", "Exon F1 [%]",
        "Tx Sens [%]",   "Tx Prec [%]",   "Tx F1 [%]"
    ]
    table = ax.table(
        cellText=cell_text,
        colLabels=col_labels,
        loc="center",
        cellLoc="center"
    )
    table.auto_set_font_size(False)
    table.set_fontsize(10)
    table.scale(1, 1.5)

def plot_gffcompare_accuracy_with_table(file_paths, labels=None):
    """
    Parse gffcompare summaries, plot Exon- & Transcript-level accuracy,
    and show a table of all metrics below the plots.

    Parameters
    ----------
    file_paths : list of str
        Paths to files containing gffcompare summaries.
    labels : list of str, optional
        Labels for each file (must match length of file_paths). If None,
        the basename of each file will be used.
    """
    # --- prepare labels
    if labels is None:
        labels = [os.path.basename(p) for p in file_paths]
    elif len(labels) != len(file_paths):
        raise ValueError("`labels` must be same length as `file_paths`")

    # --- regexes to grab Sensitivity | Precision
    exon_pat = re.compile(r'Exon level:\s*([\d\.]+)\s*\|\s*([\d\.]+)')
    tx_pat   = re.compile(r'Transcript level:\s*([\d\.]+)\s*\|\s*([\d\.]+)')

    # --- collect metrics
    exon_sens, exon_prec = [], []
    tx_sens,   tx_prec   = [], []

    for path in file_paths:
        text = open(path).read()
        m_exon = exon_pat.search(text)
        m_tx   = tx_pat.search(text)
        if not (m_exon and m_tx):
            raise RuntimeError(f"Missing Exon/Transcript lines in {path!r}")
        exon_sens.append(float(m_exon.group(1)))
        exon_prec.append(float(m_exon.group(2)))
        tx_sens.append(float(m_tx.group(1)))
        tx_prec.append(float(m_tx.group(2)))

    # --- compute F1 scores
    def f1(s, p):
        return 2 * s * p / (s + p) if (s + p) > 0 else 0.0

    exon_f1 = [f1(s, p) for s, p in zip(exon_sens, exon_prec)]
    tx_f1   = [f1(s, p) for s, p in zip(tx_sens, tx_prec)]
    # --- set up a 2x1 grid: top row for 2 plots, bottom row for the table
    fig = plt.figure(figsize=(12, 8))
    fig.suptitle("Test Species Accuracy", fontweight="bold", size=15)
    gs  = GridSpec(2, 2, height_ratios=[3, 1], hspace=0.4)

    ax1 = fig.add_subplot(gs[0, 0])
    ax2 = fig.add_subplot(gs[0, 1])
    ax_table = fig.add_subplot(gs[1, :])

    label_size=12
    title_size=14
    # Exon-level scatter
    for prec, sens, lbl in zip(exon_prec, exon_sens, labels):
        ax1.scatter(prec, sens, label=lbl, marker='x', lw=2)
    ax1.set_xlabel("Precision [%]", size=label_size)
    ax1.set_ylabel("Sensitivity [%]", size=label_size)
    ax1.set_title("Exon", fontweight="bold", size=title_size)
    ax1.set_ylim((0,100))
    ax1.set_xlim((0,100))

    # Transcript-level scatter
    for prec, sens, lbl in zip(tx_prec, tx_sens, labels):
        ax2.scatter(prec, sens, label=lbl, marker='x', lw=2)
    ax2.set_xlabel("Precision [%]", size=label_size)
    ax2.set_ylabel("Sensitivity [%]", size=label_size)
    ax2.set_title("Transcript", fontweight="bold", size=title_size)
    ax2.set_ylim((0,100))
    ax2.set_xlim((0,100))
    ax2.legend(loc="upper left")

    # now draw the table via our helper
    draw_metrics_table(
        ax_table,
        labels,
        exon_sens, exon_prec, exon_f1,
        tx_sens, tx_prec, tx_f1
    )
    plt.savefig("test_acc.png", dpi=200)
    plt.show()