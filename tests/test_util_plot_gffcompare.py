from textwrap import dedent

import matplotlib.pyplot as plt
import pytest

from tiberius_dataprep import util


_STATS = dedent(
    """\
    #-----------------| Sensitivity | Precision  |
            Base level:    88.0     |    87.3    |
            Exon level:    {exon_s}     |    {exon_p}    |
          Intron level:    74.2     |    70.3    |
    Intron chain level:    33.0     |    32.8    |
      Transcript level:    {tx_s}     |    {tx_p}    |
           Locus level:    41.7     |    36.4    |
    """
)


def _write_stats(dirpath, name, exon_s, exon_p, tx_s, tx_p):
    p = dirpath / f"{name}.stats"
    p.write_text(_STATS.format(exon_s=exon_s, exon_p=exon_p, tx_s=tx_s, tx_p=tx_p))
    return str(p)


def test_plot_gffcompare_accuracy_plots_one_point_per_file(tmp_path):
    p1 = _write_stats(tmp_path, "run1", 70.0, 65.0, 40.0, 35.0)
    p2 = _write_stats(tmp_path, "run2", 80.0, 75.0, 50.0, 45.0)

    plt.close("all")
    util.plot_gffcompare_accuracy([p1, p2], labels=["A", "B"])

    fig = plt.gcf()
    exon_ax, tx_ax = fig.axes

    exon_offsets = exon_ax.collections[0].get_offsets()
    assert exon_offsets.tolist() == [[65.0, 70.0], [75.0, 80.0]]

    tx_offsets = tx_ax.collections[0].get_offsets()
    assert tx_offsets.tolist() == [[35.0, 40.0], [45.0, 50.0]]
    plt.close("all")


def test_plot_gffcompare_accuracy_raises_on_label_mismatch(tmp_path):
    p1 = _write_stats(tmp_path, "run1", 70.0, 65.0, 40.0, 35.0)
    with pytest.raises(ValueError):
        util.plot_gffcompare_accuracy([p1], labels=["A", "B"])


def test_plot_gffcompare_accuracy_raises_on_missing_lines(tmp_path):
    bad = tmp_path / "broken.stats"
    bad.write_text("nothing useful here\n")
    with pytest.raises(RuntimeError):
        util.plot_gffcompare_accuracy([str(bad)])
