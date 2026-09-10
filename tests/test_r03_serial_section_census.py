from scripts.r03_serial_section_census import _verdicts


def test_only_in_plane_is_measurable():
    verdicts = _verdicts()
    assert verdicts["in_plane_core_masked"]["verdict"] == "MEASURABLE"
    for task in ("adjacent_plane_prediction", "held_out_middle_plane",
                 "stack_edge_extrapolation"):
        assert verdicts[task]["verdict"] == "NOT_MEASURABLE"
        assert "spacing" in verdicts[task]["reason"]
