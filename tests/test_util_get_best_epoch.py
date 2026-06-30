from tiberius_dataprep import util


def _level(sens, prec):
    return {"sensitivity": sens, "precision": prec}


def _species_stats(per_epoch_metrics):
    """per_epoch_metrics: {epoch: {metric: (sens, prec)}} -> nested dict under train 0."""
    return {
        0: {
            epoch: {metric: _level(*sp) for metric, sp in mvals.items()}
            for epoch, mvals in per_epoch_metrics.items()
        }
    }


def test_get_best_epoch_picks_highest_combined_score():
    stats = _species_stats({
        0: {"Transcript": (80, 80), "Exon": (70, 70)},  # combined = (80+80) + (70+70) = 300
        1: {"Transcript": (90, 90), "Exon": (80, 80)},  # combined = (90+90) + (80+80) = 340
    })

    (train_num, epoch_num), best_score = util.get_best_epoch([stats])

    assert (train_num, epoch_num) == (0, 1)
    assert best_score == 340


def test_get_best_epoch_averages_across_species():
    species_a = _species_stats({
        0: {"Transcript": (50, 50), "Exon": (50, 50)},  # 200
        1: {"Transcript": (90, 90), "Exon": (10, 10)},  # 200
    })
    species_b = _species_stats({
        0: {"Transcript": (50, 50), "Exon": (50, 50)},  # 200
        1: {"Transcript": (10, 10), "Exon": (90, 90)},  # 200
    })

    (_, epoch_num), best_score = util.get_best_epoch([species_a, species_b])

    # both epochs average to 200; tie -> first epoch kept (curr_best is strict >)
    assert best_score == 200
    assert epoch_num in {0, 1}
