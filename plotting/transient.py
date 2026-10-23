"""Has the growth rate settled, or is the run still in its transient?

This module exists because a measurement in this study was reported without it and was wrong by about four per cent.

A perturbed calculation does not grow at the modal rate from the first instant. The imposed shape is a single Fourier mode of the interface, but the velocity field is not the corresponding eigenvector, so the initial state is a combination of the growing mode and a set of decaying ones. Only after the decaying part has died away does the amplitude grow at the eigenvalue. Fitting an exponential before that point measures a mixture, and because the decaying part reduces the apparent rate, the fit is biased low.

The bias is not obvious in the amplitude curve, which looks like a clean exponential on a logarithmic axis long before it is one. It is obvious in the local growth rate, which rises towards a plateau; if the curve is still climbing at the end of the record, the run was stopped too early and no rate from it should be quoted. That is the check performed here, and it is the check that was missing.

The transient is modelled as a single decaying exponential,

    sigma(t) = sigma_inf - A exp(-beta t),

which is the leading behaviour once the slowest decaying mode dominates the remainder. Three quantities follow. The asymptote sigma_inf estimates what the rate would settle to; the residual A exp(-beta t_end) states how much transient is still present at the end of the record, which is the bias in the reported rate; and the time at which that residual falls below a chosen tolerance gives the run length actually required.

The model is deliberately simple and should not be trusted far beyond the data. Extrapolating an asymptote from a curve that is still steeply rising is an extrapolation, and its purpose here is diagnostic rather than quantitative: it answers whether the run is long enough, and if not roughly how much longer it must be. The rate that goes into the thesis should come from a run whose local growth rate has visibly flattened, not from this extrapolation.
"""

# ORIENTATION
#
# If you read one module in this project to understand why measurements go wrong, read this one. It exists because a rate was reported from a run that had not finished settling, and was wrong by four per cent.


from __future__ import annotations

import numpy as np

#: Fraction of the asymptote below which the remaining transient is deemed negligible. One part in a thousand is a tenth of the precision to which a growth rate is usually quoted, so a run meeting it is not limited by its length.
DEFAULT_TOLERANCE = 1.0e-3


#: Amplitude, as a fraction of the column radius, beyond which the linear description is not claimed to hold. The proportionality test in this study holds to about a tenth of the radius; past that the second-order term is no longer negligible, harmonics appear, and the mean radius drifts.
LINEAR_BOUND = 0.1


#: Fraction of the linear part of a record discarded as startup transient before a growth rate is fitted. This is the one default shared by every script that fits a rate (fit_records.py, plot_results.py), so that a rate quoted in a table and the same rate drawn in a figure come from the same window. Before it was shared the scripts used 0.5, 0.45, a one-sided half of the whole record and a hard-coded t = 4.5, and the same record gave four different numbers. Half the linear part is kept because the fit_records documentation and the archived measurement both describe discarding the first half.
DEFAULT_KEEP = 0.5


def linear_end(times, amplitudes, r0=1.0, bound=LINEAR_BOUND):
    """Last instant at which the amplitude is still small enough to be linear.

    This is the upper limit of any window a growth rate may be fitted over, and omitting it is how a correct diagnosis of one failure produces another. A run continued past this point does not merely stop being informative: the amplitude saturates and then oscillates as the column breaks and the fragments relax, so the local growth rate first exceeds the eigenvalue and then swings wildly. A window chosen as a fixed fraction of the record then straddles the break-up, and the fit returns a number that belongs to neither regime.

    Returns the time itself and whether the bound was reached at all. A record that never reaches it is entirely linear, which is the ordinary case for a short run and for every calculation in the stable band.
    """
    times = np.asarray(times, dtype=float)
    amplitudes = np.asarray(amplitudes, dtype=float)
    limit = bound * r0

    beyond = np.where(np.isfinite(amplitudes) & (amplitudes > limit))[0]
    if beyond.size == 0:
        return float(times[-1]), False
    first = int(beyond[0])
    if first == 0:
        # The very first sample already exceeds the bound, so no part of the record is linear. The caller must be told rather than handed the initial time as though it were a window.
        return float(times[0]), True
    return float(times[first - 1]), True


def select_window(times, amplitudes, r0=1.0, bound=LINEAR_BOUND,
                  keep=DEFAULT_KEEP, width=0.8):
    """Fit window bounded on both sides: after the transient, before nonlinearity.

    The upper bound comes from `linear_end`. The lower bound is placed a fixed fraction of the way through the *linear* part of the record rather than of the whole record, which is the distinction that matters: a calculation carried three times further than linearity permits has a fixed fraction of its record sitting deep in the nonlinear regime, and a window defined against the record length lands there.

    The window is returned together with a flag saying whether it is wide enough to fit over. A narrow window is not a failure of this function; it means the perturbation amplitude was too large for the transient to decay before nonlinearity arrived, and the remedy is a smaller amplitude rather than a longer run.
    """
    times = np.asarray(times, dtype=float)
    t_hi, reached = linear_end(times, amplitudes, r0=r0, bound=bound)
    t_start = float(times[0])
    t_lo = t_start + keep * (t_hi - t_start)
    wide_enough = (t_hi - t_lo) >= 2.0 * width
    return t_lo, t_hi, reached, wide_enough


def choose_window(times, amplitudes, r0=1.0, bound=LINEAR_BOUND,
                  keep=DEFAULT_KEEP, width=0.8, t_min=None, min_samples=3):
    """The fit window every script uses, with an explicit fallback and a warning.

    This wraps `select_window` so that the scripts share not only the rule but also what happens when the rule fails, which is where they used to diverge silently. An explicit `t_min` replaces the lower bound only; the upper bound stays at the end of the linear part, since no explicit choice of start time makes a fit across break-up meaningful.

    Returns (t_lo, t_hi, nonlinear_reached, warning). The warning is an empty string when the window is sound; otherwise it says what was wrong and which fallback was taken, and the caller is expected to print it beside the rate rather than discard it, because a rate from a fallback window is a weaker measurement and the reader must be told so.
    """
    times = np.asarray(times, dtype=float)
    amplitudes = np.asarray(amplitudes, dtype=float)
    t_lo, t_hi, reached, wide = select_window(times, amplitudes, r0=r0,
                                              bound=bound, keep=keep,
                                              width=width)
    t_start, t_end = float(times[0]), float(times[-1])
    if t_min is not None:
        t_lo = float(t_min)

    if reached and t_hi <= t_start:
        # No sample is linear, so there is no two-sided window at all. The older one-sided rule over the whole record is used so that a number is still produced for inspection, and the warning says plainly that it is not a linear growth rate.
        t_lo = t_min if t_min is not None else t_start + keep * (t_end - t_start)
        t_hi = t_end
        return t_lo, t_hi, reached, (
            f"no linear part (amplitude > {bound:g} r0 from the first sample); "
            f"fell back to the one-sided window t >= {t_lo:.2f} over the whole "
            "record, which is NOT a linear growth rate")

    usable = ((times >= t_lo) & (times <= t_hi) & np.isfinite(amplitudes)
              & (amplitudes > 0.0))
    if usable.sum() < min_samples:
        # Too few samples between the bounds to fit a line with any confidence. The whole linear part is used instead, which includes transient and therefore biases the rate low, and the warning says so.
        t_lo = t_start
        return t_lo, t_hi, reached, (
            f"fewer than {min_samples} samples in the window; fell back to "
            f"the whole linear part {t_start:.2f} <= t <= {t_hi:.2f}, which "
            "includes the startup transient and biases the rate low")
    if not wide:
        return t_lo, t_hi, reached, (
            f"window {t_lo:.2f} <= t <= {t_hi:.2f} is narrower than "
            f"{2.0 * width:g}; the transient may not have decayed before "
            "nonlinearity, so reduce eps rather than lengthen the run")
    return t_lo, t_hi, reached, ""


def local_growth_rate(times, amplitudes, width, n_windows=8):
    """Growth rate fitted over sliding windows of the record.

    Each window gives one estimate, obtained by fitting a straight line to the logarithm of the amplitude. Plotted against the window centre these estimates form the curve whose flattening signals that the transient has gone.
    """
    times = np.asarray(times, dtype=float)
    amplitudes = np.asarray(amplitudes, dtype=float)

    lo = times.min() + 0.5 * width
    hi = times.max() - 0.5 * width
    if hi <= lo:
        return np.array([]), np.array([])

    centres, rates = [], []
    for centre in np.linspace(lo, hi, n_windows):
        inside = (times >= centre - 0.5 * width) & (times <= centre + 0.5 * width)
        inside &= np.isfinite(amplitudes) & (amplitudes > 0.0)
        if inside.sum() >= 4:
            # A straight line through log a has slope sigma, since a = a0 exp(sigma t) implies log a = log a0 + sigma t.
            slope = np.polyfit(times[inside], np.log(amplitudes[inside]), 1)[0]
            centres.append(centre)
            rates.append(slope)
    return np.asarray(centres), np.asarray(rates)


def fit_approach(centres, rates):
    """Fit sigma(t) = sigma_inf - A exp(-beta t) to the local rate.

    Returns the three parameters, or None if the fit does not converge, which happens when the record is too short or too noisy to constrain three unknowns. Returning None rather than an unconverged answer matters, because the caller uses this to decide whether a measurement may be quoted.
    """
    from scipy.optimize import curve_fit

    centres = np.asarray(centres, dtype=float)
    rates = np.asarray(rates, dtype=float)
    if centres.size < 4:
        return None

    # The fit is made in time measured from the first window centre rather than in absolute time. In absolute time the amplitude parameter is the transient extrapolated back to t = 0, which for a window starting at t = 5 with beta = 1 is e^5 times larger than anything in the data and strongly correlated with beta, so the least-squares problem is badly conditioned and the starting guess below means nothing. Shifted, the amplitude is the transient at the first centre, of the order of the observed rise.
    t_ref = float(centres[0])
    shifted = centres - t_ref

    def model(s, sigma_inf, amp, beta):
        return sigma_inf - amp * np.exp(-beta * s)

    # The starting guess matters, and it must allow for either direction of approach. A rate rising to its plateau has its asymptote just above the last estimate and a positive amplitude; a rate falling to it from above has its asymptote just below and a negative amplitude. The earlier guess assumed the rising case only, so an approach from above never converged and was reported as unfittable. The amplitude at the first centre follows from the model at s = 0, sigma = sigma_inf - A, and the decay rate is of order the reciprocal of the record length.
    # ndarray.ptp was removed in NumPy 2; the free function is the portable spelling.
    span = float(np.ptp(rates))
    direction = 1.0 if rates[-1] >= rates[0] else -1.0
    sigma_guess = float(rates[-1]) + direction * 0.05 * max(span, 1e-12)
    guess = [sigma_guess,
             sigma_guess - float(rates[0]),
             1.0 / max(float(np.ptp(centres)), 1.0)]
    try:
        params, _ = curve_fit(model, shifted, rates, p0=guess, maxfev=20000)
    except (RuntimeError, ValueError):
        return None

    sigma_inf, amp_ref, beta = (float(v) for v in params)
    if beta <= 0.0 or not np.isfinite(sigma_inf):
        # A non-positive decay rate means the model has fitted a growing transient, which is not the situation this diagnostic describes.
        return None
    # Converted back to the absolute-time amplitude A of sigma_inf - A exp(-beta t), so that the returned parameters keep the meaning stated in the module docstring and residual_transient and required_time can still be evaluated at absolute times.
    amp = amp_ref * np.exp(beta * t_ref)
    if not np.isfinite(amp):
        return None
    return sigma_inf, float(amp), beta


def residual_transient(amp, beta, t_end):
    """Transient still present at the end of the record, in absolute units."""
    return float(amp * np.exp(-beta * t_end))


def required_time(amp, beta, sigma_inf, tolerance=DEFAULT_TOLERANCE):
    """Time at which the remaining transient falls below a tolerance.

    Obtained by inverting |A| exp(-beta t) = tolerance . |sigma_inf|. The magnitude of A is used because a negative A is a rate approaching its plateau from above, which is as much a transient as an approach from below and decays at the same rate; returning NaN for it, as this once did, hid the run length such a case needs.
    """
    target = tolerance * abs(sigma_inf)
    if target <= 0.0 or amp == 0.0 or not np.isfinite(amp):
        return float("nan")
    return float(np.log(abs(amp) / target) / beta)


def assess(times, amplitudes, width=0.8, tolerance=DEFAULT_TOLERANCE,
           sigma_reference=None, r0=1.0, bound=LINEAR_BOUND):
    """Full diagnosis of whether a record is long enough to quote a rate.

    The returned dictionary carries `settled`, which is the answer to the question the caller is really asking. Everything else supports it: the local rate curve for plotting, the fitted asymptote, the bias remaining in the record, and the run length that would remove it.
    """
    times = np.asarray(times, dtype=float)
    amplitudes = np.asarray(amplitudes, dtype=float)

    # The record is cut back to its linear part before anything is fitted. Assessing the settling of a transient on a record that includes break-up fits the model to the wrong thing entirely, and the earlier version of this function did exactly that: it returned an asymptote above the inviscid limit, which is impossible, because the nonlinear rise dominated the fit.
    t_hi, reached = linear_end(times, amplitudes, r0=r0, bound=bound)
    linear = times <= t_hi
    times_lin = times[linear]
    amps_lin = amplitudes[linear]

    centres, rates = local_growth_rate(times_lin, amps_lin, width)

    result = {
        "centres": centres,
        "rates": rates,
        "settled": False,
        "sigma_inf": float("nan"),
        "beta": float("nan"),
        "residual": float("nan"),
        "relative_bias": float("nan"),
        "t_required": float("nan"),
        "nonlinear_reached": reached,
        "t_linear_end": t_hi,
        "message": "",
    }
    if reached and t_hi <= float(times[0]):
        result["message"] = (
            "THE PERTURBATION IS TOO LARGE. The amplitude exceeds "
            f"{bound:g} r0 from the first sample, so no part of this record "
            "is linear and no growth rate can be fitted. Reduce eps.")
        return result
    if centres.size < 4:
        if reached:
            result["message"] = (
                "THE PERTURBATION IS TOO LARGE. The amplitude leaves the "
                f"linear range at t = {t_hi:.2f}, which is too early for the "
                "transient to have decayed, so there is no window in which a "
                "rate may be fitted. Reduce eps rather than lengthening the "
                "run.")
        else:
            result["message"] = ("record too short to assess; fewer than four "
                                 "sliding windows fit inside it")
        return result

    # The final rise, expressed against the spread of the whole curve, is a model-free indicator. A curve that has flattened shows a last increment far smaller than its total rise.
    total_rise = float(rates.max() - rates.min())
    last_step = float(rates[-1] - rates[-2])
    result["last_increment"] = last_step
    result["total_rise"] = total_rise

    fit = fit_approach(centres, rates)
    if fit is None:
        result["message"] = ("the approach to a plateau could not be fitted; "
                             "inspect the local growth rate by eye")
        return result

    sigma_inf, amp, beta = fit
    # The end of the linear part, not the end of the record, since the rate is only ever quoted from inside the linear part.
    t_end = float(times_lin.max())
    residual = residual_transient(amp, beta, t_end)

    result.update({
        "sigma_inf": sigma_inf,
        "beta": beta,
        "residual": residual,
        "relative_bias": residual / abs(sigma_inf) if sigma_inf else np.nan,
        "t_required": required_time(amp, beta, sigma_inf, tolerance),
        # The magnitude of the residual is compared, not its signed value. A rate approaching its plateau from above has a negative fitted amplitude, hence a negative residual, which a signed comparison would always pass however large the remaining transient. The relative bias above keeps its sign, since the direction of the bias is worth reporting.
        "settled": abs(residual) / abs(sigma_inf) < tolerance if sigma_inf else False,
    })

    if result["settled"]:
        note = ""
        if reached:
            note = (f" The record continues past t = {t_hi:.2f}, where the "
                    "amplitude leaves the linear range; that part is excluded.")
        result["message"] = (
            f"the local growth rate has settled: the remaining transient is "
            f"{100 * result['relative_bias']:.2f} per cent of the asymptote, "
            f"below the {100 * tolerance:.2f} per cent tolerance." + note)
    elif reached:
        # Both bounds bite at once, which is the situation a longer run cannot fix.
        result["message"] = (
            f"NO USABLE WINDOW. The transient still contributes about "
            f"{100 * result['relative_bias']:.2f} per cent at t = {t_hi:.2f}, "
            f"where the amplitude already leaves the linear range. Lengthening "
            f"the run will not help, because the amplitude grows past "
            f"{bound:g} r0 before the transient decays. Reduce eps: the "
            f"transient needs until t = {result['t_required']:.2f}, so eps "
            f"must be small enough that the amplitude is still below "
            f"{bound:g} r0 then.")
    else:
        result["message"] = (
            f"THE RUN IS TOO SHORT. The local growth rate is still rising, and "
            f"about {100 * result['relative_bias']:.2f} per cent of transient "
            f"remains at t = {t_end:.2f}, which biases the fitted rate low by "
            f"roughly that amount. Advance the calculation to t = "
            f"{result['t_required']:.2f} before quoting a rate.")

    if sigma_reference is not None and np.isfinite(sigma_inf):
        result["reference_gap"] = (sigma_inf - sigma_reference) / sigma_reference

    return result


def n_periods_required(t_required, sigma_reference):
    """Run length in e-folding times, the unit the deck generator expects.

    `make_inputs` measures the stop time in e-folding times of the reference rate, so a required time is converted here rather than in the caller, where the convention would have to be restated.
    """
    if not np.isfinite(t_required) or sigma_reference <= 0.0:
        return float("nan")
    return float(t_required * sigma_reference)


def format_report(result, sigma_reference=None):
    """The diagnosis as plain text, for printing beneath a measurement."""
    lines = ["-" * 66, "transient convergence", "-" * 66]
    if np.isfinite(result.get("sigma_inf", np.nan)):
        lines.append(f"{'extrapolated asymptote':<28} "
                     f"{result['sigma_inf']:.6f}")
        lines.append(f"{'transient decay rate beta':<28} "
                     f"{result['beta']:.4f}")
        lines.append(f"{'transient left in the record':<28} "
                     f"{100 * result['relative_bias']:.2f} per cent")
        lines.append(f"{'time needed for convergence':<28} "
                     f"{result['t_required']:.2f}")
        if sigma_reference:
            lines.append(f"{'as e-folding times':<28} "
                         f"{n_periods_required(result['t_required'], sigma_reference):.2f}")
    lines.append("")
    lines.append(result["message"])
    lines.append("-" * 66)
    return "\n".join(lines)
