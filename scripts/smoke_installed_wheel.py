"""Install a wheel into a temporary target and probe it outside the checkout."""

import argparse
import os
import subprocess
import sys
import tempfile
from pathlib import Path

PROBE = r'''
import importlib.util
from importlib import resources
import numpy as np
from pathlib import Path

import stiffgwpy_fast
from stiffgwpy_fast import LCDM_SG

required = {
    "stiffgwpy_fast": ("th.dat", "fd_table.npz"),
}
if importlib.util.find_spec("cobaya") is not None:
    # 安装可选 Cobaya 依赖时，同时核验其插件和似然数据资源。
    required.update({
    "stiffgwpy_fast.cobaya": ("stiffGW.yaml",),
    "stiffgwpy_fast.cobaya.likelihoods.LIGO_SGWB": ("C_O1_O2_O3.dat", "LVK_SGWB_CC.yaml"),
    "stiffgwpy_fast.cobaya.likelihoods.PTA": ("IPTA.yaml", "NANOGrav.yaml"),
    "stiffgwpy_fast.cobaya.likelihoods.PTA.EPTAdr2": ("EPTA_dr2new_mock.dat", "freqs_dr2new.txt"),
    "stiffgwpy_fast.cobaya.likelihoods.PTA.NANOGrav15yr": (
        "density_mock.npy", "freqs.npy", "log10rhogrid.npy"
    ),
    })
for package, names in required.items():
    root = resources.files(package)
    for name in names:
        item = root.joinpath(name)
        assert item.is_file(), f"missing package resource: {package}/{name}"
        with resources.as_file(item) as path:
            assert Path(path).is_file()

with resources.as_file(resources.files("stiffgwpy_fast").joinpath("th.dat")) as path:
    table = np.loadtxt(path)
assert table.ndim == 2 and table.shape[1] >= 7
model = LCDM_SG(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2)
assert model.derived_param["N_inf"] is not None
print("installed-wheel smoke: OK", stiffgwpy_fast.__version__)
'''


def main(argv=None):
    parser = argparse.ArgumentParser()
    parser.add_argument("wheel", type=Path)
    args = parser.parse_args(argv)
    wheel = args.wheel.resolve()
    if wheel.suffix != ".whl" or not wheel.is_file():
        raise SystemExit(f"not a wheel file: {wheel}")

    with tempfile.TemporaryDirectory(prefix="stiffgwpy_fast-wheel-smoke-") as tmp:
        root = Path(tmp)
        target = root / "site"
        cwd = root / "outside-checkout"
        target.mkdir()
        cwd.mkdir()
        subprocess.check_call([
            sys.executable, "-m", "pip", "install", "--disable-pip-version-check",
            "--no-deps", "--target", os.fspath(target), os.fspath(wheel),
        ])
        env = os.environ.copy()
        env["PYTHONPATH"] = os.fspath(target)
        subprocess.check_call([sys.executable, "-c", PROBE], cwd=cwd, env=env)


if __name__ == "__main__":
    main()
