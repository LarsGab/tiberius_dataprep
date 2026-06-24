from tiberius_dataprep import util


def test_read_training_stats_empty_species_list_returns_empty_dict(tmp_path):
    assert util.read_training_stats(str(tmp_path), species_list=[]) == {}
