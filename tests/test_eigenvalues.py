import hashlib, os, pathlib
from fractions import Fraction
import pytest
from bolzakit import spectral
from bolzakit.spectral import eigen as E, zeta
from bolzakit.errors import DataUnavailable
from bolzakit.provenance import info


def test_fem_clusters():
    c = E.eigenvalues("fem")
    assert [m for _, m in c[:5]] == [3, 4, 2, 4, 3]
    assert abs(c[0][0] - 3.8388872588421995) / 3.8388872588421995 < 2e-4


def test_fem_beyond_bundle():
    with pytest.raises(DataUnavailable, match=r"FEM bundle has \d+ clusters"):
        E.eigenvalues("fem", n=100)


def test_zeta_zero():
    assert zeta.zeta_zero() == Fraction(-4, 3)
    assert info(zeta.zeta_zero).label == "PROVED"


def test_su_corrupt_cache(tmp_path, monkeypatch):
    monkeypatch.setenv("BOLZAKIT_CACHE", str(tmp_path))
    bad = tmp_path / "eig-bolza-refined0-1000.txt"
    bad.write_text("3.0\n")
    monkeypatch.setattr(E, "_download", lambda dest: dest.write_text("garbage\n"))
    with pytest.raises(DataUnavailable, match="checksum"):
        E.su2013()
    assert not bad.exists()


def test_su_offline_message(tmp_path, monkeypatch):
    monkeypatch.setenv("BOLZAKIT_CACHE", str(tmp_path))
    def boom(dest):
        raise OSError("no network")
    monkeypatch.setattr(E, "_download", boom)
    with pytest.raises(DataUnavailable, match="arxiv.org"):
        E.su2013()


@pytest.mark.network
def test_su_real_download_and_agreement_with_fem():
    vals = E.su2013()
    assert len(vals) == 1000
    su = E.eigenvalues("su2013")
    for (vf, mf), (vs, ms) in zip(E.eigenvalues("fem")[:5], su[:5]):
        assert mf == ms and abs(vf - vs) / vs < 2e-4


def test_spectral_eigenvalues_is_callable():
    assert callable(spectral.eigenvalues) and spectral.eigenvalues("fem")[0][1] == 3


def test_su_explicit_bad_path_kept(tmp_path):
    bad = tmp_path / "mine.txt"
    bad.write_text("3.0\n")
    with pytest.raises(DataUnavailable, match="checksum"):
        E.su2013(path=str(bad))
    assert bad.exists()


def test_su_explicit_missing_path(tmp_path):
    with pytest.raises(DataUnavailable):
        E.su2013(path=str(tmp_path / "nope.txt"))


def test_eigenvalues_n_zero():
    with pytest.raises(ValueError):
        E.eigenvalues("fem", n=0)


def test_su_eigenvalues_alias():
    assert callable(E.su_eigenvalues)
