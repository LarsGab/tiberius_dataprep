from textwrap import dedent

from tiberius_dataprep import util


_STATS_TEMPLATE = dedent(
    """\
    #-----------------| Sensitivity | Precision  |
            Base level:    {base_s}     |    {base_p}    |
            Exon level:    {exon_s}     |    {exon_p}    |
          Intron level:    74.2     |    70.3    |
    Intron chain level:    33.0     |    32.8    |
      Transcript level:    {tx_s}     |    {tx_p}    |
           Locus level:    41.7     |    36.4    |
    """
)


def _write_stats(dirpath, train_num, epoch_num, species, **levels):
    fname = dirpath / f"{train_num}.epoch_{epoch_num:02d}.{species}.stats"
    fname.write_text(_STATS_TEMPLATE.format(**levels))


def _full_levels(base_s, base_p, exon_s, exon_p, tx_s, tx_p):
    return dict(
        base_s=base_s, base_p=base_p,
        exon_s=exon_s, exon_p=exon_p,
        tx_s=tx_s, tx_p=tx_p,
    )


def test_read_training_stats_parses_two_epochs(tmp_path):
    _write_stats(tmp_path, 0, 1, "athal", **_full_levels(88.0, 87.3, 68.8, 64.9, 38.6, 36.4))
    _write_stats(tmp_path, 0, 2, "athal", **_full_levels(90.0, 89.0, 70.0, 66.0, 42.0, 40.0))

    out = util.read_training_stats(str(tmp_path), species_list=["athal"])

    assert set(out.keys()) == {"athal"}
    assert set(out["athal"].keys()) == {0}
    assert set(out["athal"][0].keys()) == {1, 2}
    assert out["athal"][0][1]["Exon"] == {"sensitivity": 68.8, "precision": 64.9}
    assert out["athal"][0][2]["Transcript"] == {"sensitivity": 42.0, "precision": 40.0}


def test_read_training_stats_ignores_other_species(tmp_path):
    _write_stats(tmp_path, 0, 1, "athal", **_full_levels(88.0, 87.3, 68.8, 64.9, 38.6, 36.4))
    _write_stats(tmp_path, 0, 1, "oryza", **_full_levels(70.0, 60.0, 50.0, 40.0, 30.0, 20.0))

    out = util.read_training_stats(str(tmp_path), species_list=["athal"])

    assert "oryza" not in out
    assert out["athal"][0][1]["Base"]["sensitivity"] == 88.0
