from stiffgwpy_fast import fast_sgwb as FS
from stiffgwpy_fast.stiff_SGWB import LCDM_SG


def test_exact_frequency_preparation_reuses_model_workspace():
    model = LCDM_SG(r=1e-2, cr=1, T_re=2e3, kappa10=1e-2)
    FS.gen_fast(model, kink_split=True)
    if not hasattr(model, 'f'):
        model.construct_f(1.0)
    freqs = model.f
    first = FS.prep_frequency_only(model, model.Nv, freqs)
    second = FS.prep_frequency_only(model, model.Nv, freqs)

    for before, after in zip(first[2:], second[2:]):
        assert before is after
