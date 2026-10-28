import subprocess
import sys
from pathlib import Path


def test_pure_pta_helper_import_does_not_require_optional_cobaya():
    root = Path(__file__).resolve().parents[1]
    code = """
import builtins

original_import = builtins.__import__

def reject_cobaya(name, *args, **kwargs):
    if name == 'cobaya' or name.startswith('cobaya.'):
        raise ModuleNotFoundError('Cobaya is intentionally unavailable')
    return original_import(name, *args, **kwargs)

builtins.__import__ = reject_cobaya
from stiffgwpy_fast.cobaya.likelihoods.PTA._spectrum import primordial_rho
assert callable(primordial_rho)
"""
    result = subprocess.run(
        [sys.executable, "-c", code],
        cwd=root,
        capture_output=True,
        text=True,
        check=False,
    )

    assert result.returncode == 0, result.stderr
