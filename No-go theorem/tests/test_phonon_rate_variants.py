"""Guards the Gamma_XY(T) model-variant selector (phonon_rates.py) added for
the room-temperature no-go plan (Sec. 4.1: "no-goに最も不利な、すなわち散逸を
最も小さく見積もるモデルを主結果に用いる" -- use the model that estimates the
SMALLEST dissipation as the main result).

The key failure mode this guards against: the "saturation" variant's ceiling
is a free construction parameter (not a fitted literature quantity, unlike
the other three variants). If chosen carelessly it can silently win the
"most conservative" selection everywhere by construction alone, which would
make the no-go campaign's main result rest on an arbitrary number instead of
a genuine physical bound. This test pins that the ceiling does NOT bind
within the declared campaign range T in [4, 300] K (SMRT plan Sec. 4.1), so
"most conservative" reflects the literature-uncertainty lower bound there,
not the saturation construction.
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))
import phonon_rates as pr

CAMPAIGN_TEMPS = [4, 10, 30, 50, 77, 150, 300]
D_PERP_GHZ = 1.683  # candidate-point strain scale, matches gate scripts


def test_saturation_ceiling_does_not_bind_in_campaign_range():
    for T in CAMPAIGN_TEMPS:
        variants = pr.k_orb_variants(T, D_PERP_GHZ)
        assert variants["saturation"] > 0.9 * variants["full_happacher"], (
            f"T={T} K: saturation ceiling is binding inside the declared "
            f"campaign range (saturation={variants['saturation']:.3e} Hz vs "
            f"full_happacher={variants['full_happacher']:.3e} Hz) -- the "
            f"ceiling is acting as an unphysical shortcut, not a stress-test "
            f"cap; raise GAMMA_SAT_CEILING_HZ."
        )


def test_most_conservative_is_literature_lower_bound_in_campaign_range():
    for T in CAMPAIGN_TEMPS:
        label, val, variants = pr.k_orb_most_conservative(T, D_PERP_GHZ)
        assert label == "conservative_lower_bound", (
            f"T={T} K: expected the literature-uncertainty lower bound to "
            f"be the most-conservative (smallest dissipation) variant in "
            f"the declared campaign range, got '{label}'."
        )
        assert val == min(variants.values())


if __name__ == "__main__":
    test_saturation_ceiling_does_not_bind_in_campaign_range()
    test_most_conservative_is_literature_lower_bound_in_campaign_range()
    print("PASS: saturation ceiling does not silently dominate the campaign-range selection.")
