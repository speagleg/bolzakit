import mpmath as mp
import pytest
from bolzakit import green


def test_closed_forms():
    c = green.constants(dps=30)
    with mp.workdps(30):                        # the 20-digit literals and 1/(8 pi) must not be rounded to 15 digits
        assert mp.almosteq(c["B_Sigma_00"]["value"], mp.mpf("0.69660196484283842959"), 1e-19)
        assert mp.almosteq(c["beta1"]["value"], mp.mpf("-0.0945730196647619"), 1e-15)
        assert mp.almosteq(c["beta2"]["value"], 1 / (8 * mp.pi), 1e-28)
        assert mp.almosteq(c["beta3"]["value"], mp.mpf("-0.16248381288685"), 1e-13)


def test_lambda3_cubic_coefficient_matches_item03():
    c = green.constants()
    coef = mp.mpf(-61) / 2 + 126 * mp.pi * c["beta3"]["value"]
    assert abs(coef - mp.mpf("-94.8177")) < 1e-3


def test_bergman_numeric_matches_closed_form():
    v = green.bergman_at_weierstrass()          # dps=50, M=24, npts=50: the item 03 driver's first setting
    assert abs(v - mp.mpf("0.69660196484283842959")) < 1e-10


def test_constants_provenance_cites_derivation():
    from bolzakit.provenance import info
    src = info(green.constants).source
    assert "review_1" in src and "docs/green.md" in src


def test_constants_match_reference_json():
    # final review M7: data/green_constants.json is a reference copy of the closed form, not a source of values
    from bolzakit import data
    ref = data.load("green_constants")
    c = green.constants(dps=30)
    assert set(ref) == set(c)
    with mp.workdps(40):
        for k in c:
            assert abs(mp.mpf(ref[k]["value"]) - c[k]["value"]) < mp.mpf("1e-25"), k
            assert ref[k]["closed_form"] == c[k]["closed_form"] and ref[k]["label"] == c[k]["label"]
