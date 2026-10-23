import numpy as np
from scipy.interpolate import CubicSpline


def primordial_rho(f_theory, log10_omega_theory, work_freqs, hubble, time_baseline):
    rho = np.zeros_like(work_freqs, dtype=np.float64)
    theory_freqs = np.asarray(f_theory, dtype=np.float64)
    theory_omega = np.asarray(log10_omega_theory, dtype=np.float64)
    data_freqs = np.asarray(work_freqs, dtype=np.float64)

    signal_band = (theory_freqs > -13.0) & (theory_freqs < -5.0)
    knots = theory_freqs[signal_band]
    values = theory_omega[signal_band]
    if knots.size == 0:
        return rho

    query = np.log10(data_freqs)
    if knots.size == 1:
        supported = query == knots[0]
        if not np.any(supported):
            return rho
        omega = values[0]
    else:
        supported = (query >= knots[0]) & (query <= knots[-1])
        if not np.any(supported):
            return rho
        spectrum = CubicSpline(knots, values, extrapolate=False)
        omega = spectrum(query[supported])

    rho[supported] = (
        np.power(10.0, omega) * hubble**2
        / (8.0 * np.pi**4 * data_freqs[supported]**5 * time_baseline)
    )
    return rho
