import matplotlib.pyplot as plt
import pytest
import scipy.stats

from tiberius_dataprep import util


def _stats(sens_tx, prec_tx, sens_exon, prec_exon):
    return {
        0: {
            "Transcript": {"sensitivity": sens_tx, "precision": prec_tx},
            "Exon": {"sensitivity": sens_exon, "precision": prec_exon},
        }
    }


def test_plot_accuracies_uses_harmonic_mean_for_exon():
    # Regression for util.py bug where exon accuracy was computed as the
    # SUM of sensitivity + precision instead of their harmonic mean.
    sens_exon, prec_exon = 90.0, 70.0
    plt.close("all")
    util.plot_accuracies([_stats(80.0, 80.0, sens_exon, prec_exon)])

    exon_ax = plt.gcf().axes[1]
    plotted = exon_ax.lines[0].get_ydata()[0]
    expected = scipy.stats.hmean([sens_exon, prec_exon])

    assert plotted == pytest.approx(expected)
    plt.close("all")


def test_plot_accuracies_uses_harmonic_mean_for_transcript():
    sens_tx, prec_tx = 80.0, 60.0
    plt.close("all")
    util.plot_accuracies([_stats(sens_tx, prec_tx, 90.0, 70.0)])

    tx_ax = plt.gcf().axes[0]
    plotted = tx_ax.lines[0].get_ydata()[0]
    expected = scipy.stats.hmean([sens_tx, prec_tx])

    assert plotted == pytest.approx(expected)
    plt.close("all")
