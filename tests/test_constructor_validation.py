import pytest
import yaml

from stiffgwpy_fast import LCDM_SN


def test_constructor_rejects_unknown_keyword_instead_of_using_default():
    with pytest.raises(TypeError, match="Omega_bh22"):
        LCDM_SN(Omega_bh22=0.01)


def test_constructor_rejects_unknown_key_in_parameter_mapping():
    params = LCDM_SN().cosmo_param
    params["Omega_bh22"] = 0.01

    with pytest.raises(TypeError, match="Omega_bh22"):
        LCDM_SN(params)


def test_constructor_reports_missing_yaml_file(tmp_path):
    missing = tmp_path / "missing.yaml"

    with pytest.raises(FileNotFoundError, match="missing.yaml"):
        LCDM_SN(missing)


def test_constructor_accepts_pathlike_yaml_mapping(tmp_path):
    expected = LCDM_SN().cosmo_param
    config_path = tmp_path / "model.yaml"
    config_path.write_text(yaml.safe_dump(expected), encoding="utf-8")

    model = LCDM_SN(config_path)

    assert model.cosmo_param == expected


def test_constructor_rejects_yaml_without_required_parameters(tmp_path):
    config_path = tmp_path / "partial.yaml"
    config_path.write_text("H0: 67.3\n", encoding="utf-8")

    with pytest.raises(ValueError, match="missing required"):
        LCDM_SN(config_path)


def test_constructor_surfaces_yaml_syntax_errors(tmp_path):
    config_path = tmp_path / "invalid.yaml"
    config_path.write_text("H0: [67.3\n", encoding="utf-8")

    with pytest.raises(yaml.YAMLError):
        LCDM_SN(config_path)
