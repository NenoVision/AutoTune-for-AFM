import numpy as np
from scipy import sparse
from scipy.sparse.linalg import spsolve
from scipy.signal import find_peaks
from numpy.lib.stride_tricks import sliding_window_view
import pywt
from typing import Optional

FloatArray = np.ndarray
Number = float


def baseline_als(y: FloatArray, lam: Number, p: Number, niter: int = 10) -> FloatArray:
    """
    Perform baseline correction using Asymmetric Least Squares (ALS) smoothing.

    Args:
        y (FloatArray): Input signal.
        lam (Number): Smoothness parameter.
        p (Number): Asymmetry parameter.
        niter (int): Number of iterations.

    Returns:
        FloatArray: Estimated baseline.
    """
    L = len(y)
    D = sparse.diags([1, -2, 1], [0, -1, -2], shape=(L, L - 2))
    w = np.ones(L)
    for _ in range(niter):
        W = sparse.spdiags(w, 0, L, L)
        Z = W + lam * D @ D.T
        z = spsolve(Z, w * y)
        w = p * (y > z) + (1 - p) * (y < z)
    return z


def remove_baseline(
    signal: FloatArray, lam: Number = 1_000_000, p: Number = 0.08
) -> FloatArray:
    """
    Remove baseline from a signal using ALS smoothing.

    Args:
        signal (FloatArray): Input signal.
        lam (Number): Smoothness parameter.
        p (Number): Asymmetry parameter.

    Returns:
        FloatArray: Detrended signal.
    """
    baseline = baseline_als(signal, lam, p)
    return signal - baseline


def normalize(signal: FloatArray) -> FloatArray:
    """
    Normalize a signal to the range [0, 1].

    Args:
        signal (FloatArray): Input signal.

    Returns:
        FloatArray: Normalized signal.
    """
    return (signal - np.min(signal)) / (np.max(signal) - np.min(signal))


def signal_variance(
    signal: FloatArray, window: int = 31, keep_length: bool = False
) -> FloatArray:
    """
    Compute the variance of a signal in a moving window.

    Args:
        signal (FloatArray): Input signal.
        window (int): Window size (must be odd).
        keep_length (bool): Pad the result to match input length.

    Returns:
        FloatArray: Variance signal.
    """
    variance = sliding_window_view(signal, window).var(axis=1)
    if keep_length:
        pad = (window - 1) // 2
        start = np.full(pad, variance[0])
        end = np.full(pad, variance[-1])
        variance = np.concatenate([start, variance, end])
    return variance


def find_peak(signal: FloatArray, prominence: Number = 0.85) -> Optional[int]:
    """
    Find a single dominant peak in a signal.

    Args:
        signal (FloatArray): Input signal.
        prominence (Number): Minimum relative prominence.

    Returns:
        int | None: Index of the found peak, or None if not found.
    """
    peaks, _ = find_peaks(signal, prominence=signal.max() * prominence)
    if len(peaks) == 0:
        return None
    return peaks[np.argmax(signal[peaks])]


def find_peak_bounds(
    signal: FloatArray,
    freqs: FloatArray,
    peak_index: int,
    threshold: Number = 0.01,
    window: int = 7,
) -> tuple[Number, Number]:
    """
    Determine the frequency range around a peak based on signal gradient variance.

    Args:
        signal (FloatArray): Input amplitude signal.
        freqs (FloatArray): Corresponding frequencies.
        peak_index (int): Index of the peak.
        threshold (Number): Variance threshold.
        window (int): Variance window size.

    Returns:
        tuple[Number, Number]: Minimum and maximum frequency around the peak.
    """
    gradient = np.gradient(signal)
    var = signal_variance(gradient, window, keep_length=True)
    norm_var = normalize(var)
    filtered = norm_var < threshold
    left = filtered[:peak_index][::-1]
    right = filtered[peak_index:]
    left_index = peak_index - 1 - np.argmax(left) if left.any() else 0
    right_index = peak_index + 1 + np.argmax(right) if right.any() else len(signal) - 1
    return freqs[left_index], freqs[right_index]


def _calculate_wavelet_coeffs(
    signal: FloatArray,
    frequency_step: float,
    wavelet: str = "mexh",
    widths: np.ndarray = np.arange(1, 100),
) -> np.ndarray:
    """
    Compute the wavelet transform coefficients.

    Args:
        signal (FloatArray): Input signal.
        frequency_step (float): Frequency sampling step.
        wavelet (str): Wavelet type.
        widths (np.ndarray): Wavelet widths.

    Returns:
        np.ndarray: CWT coefficients.
    """
    return pywt.cwt(signal, widths * (100 / frequency_step), wavelet)[0]


def _find_local_maxima(coeffs: np.ndarray, window_size: int = 200) -> np.ndarray:
    """
    Identify local maxima in the CWT coefficient matrix.

    Args:
        coeffs (np.ndarray): Wavelet coefficients.
        window_size (int): Size of the detection window.

    Returns:
        np.ndarray: Matrix with normalized local maxima.
    """
    local_maxima = np.zeros_like(coeffs)
    half_window = window_size // 2
    for col in range(half_window, coeffs.shape[1] - half_window):
        window_indices = np.argmax(
            coeffs[:, col - half_window : col + half_window], axis=1
        )
        local_maxima[
            np.arange(coeffs.shape[0]), col - half_window + window_indices
        ] += 1
    return _normalize_matrix(local_maxima, window_size)


def _normalize_matrix(matrix: np.ndarray, norm_factor: Number) -> np.ndarray:
    """
    Normalize a matrix by a scalar factor.

    Args:
        matrix (np.ndarray): Input matrix.
        norm_factor (Number): Normalization factor.

    Returns:
        np.ndarray: Normalized matrix.
    """
    return matrix / norm_factor


def _find_local_lines(
    local_maxima: np.ndarray, value_threshold: float = 0.5, max_index_offset: int = 100
) -> list[Number]:
    """
    Trace lines through normalized local maxima matrix.

    Args:
        local_maxima (np.ndarray): Normalized maxima matrix.
        value_threshold (float): Minimum value to consider a peak.
        max_index_offset (int): Maximum deviation allowed.

    Returns:
        list[Number]: Detected peak indices.
    """

    def _calc_window_size(
        x,
        a=-2.0874416652577305e-05,
        b=0.0011971137077197358,
        c=0.4901163045880604,
        d=85.0,
    ):
        return int(a * x**3 + b * x**2 + c * x + d)

    return_indices = []
    for wave_pos_idx in range(local_maxima.shape[1]):
        start_idx = wave_pos_idx
        width_idx = 0
        pos_array = []
        while (
            width_idx < local_maxima.shape[0]
            and local_maxima[width_idx, wave_pos_idx] >= value_threshold
            and np.abs(start_idx - wave_pos_idx) < max_index_offset
        ):
            width_idx += 1
            left = max(0, wave_pos_idx - _calc_window_size(width_idx))
            right = min(
                local_maxima.shape[1], wave_pos_idx + _calc_window_size(width_idx)
            )
            window_vals = local_maxima[width_idx, left:right]
            max_idx = np.argmax(window_vals)
            wave_pos_idx = left + max_idx
            pos_array.append(wave_pos_idx)
            if width_idx == local_maxima.shape[0] - 1:
                return_indices.append(np.mean(pos_array))
                break
    return return_indices


def _filter_peaks(
    local_lines: list[Number],
    signal: FloatArray,
    min_height: float = 0.02,
    window_size: int = 80,
) -> list[Number]:
    """
    Filter out low-amplitude or ambiguous peaks.

    Args:
        local_lines (list[Number]): Detected peak indices.
        signal (FloatArray): Original signal.
        min_height (float): Minimum required height difference.
        window_size (int): Range around peak for evaluation.

    Returns:
        list[Number]: Validated peak indices.
    """
    filtered = []
    length = len(signal)
    for idx in local_lines:
        start = max(0, int(idx - window_size))
        end = min(length, int(idx + window_size))
        window = signal[start:end]
        peak_idx = np.argmax(window)
        if (
            window[peak_idx] - window[0] > min_height
            and window[peak_idx] - window[-1] > min_height
        ):
            filtered.append(peak_idx + start)
    return filtered


def detect_peaks(
    signal_freqs: FloatArray,
    signal_amps: FloatArray,
    wavelet: str = "mexh",
    widths: np.ndarray = np.arange(1, 100),
) -> dict[float, tuple[float, float]]:
    """
    Detect and localize peaks in a frequency signal using wavelet analysis.

    Args:
        signal_freqs (FloatArray): Array of frequency values.
        signal_amps (FloatArray): Array of amplitude values.
        wavelet (str): Wavelet function to use.
        widths (np.ndarray): Wavelet widths.

    Returns:
        dict[float, tuple[float, float]]: Mapping from peak frequency to (freq_min, freq_max).
    """
    frequency_step = signal_freqs[1] - signal_freqs[0]
    detrended = remove_baseline(signal_amps)
    coeffs = _calculate_wavelet_coeffs(detrended, frequency_step, wavelet, widths)
    maxima = _find_local_maxima(coeffs)
    lines = _find_local_lines(maxima)
    filtered = _filter_peaks(lines, detrended)
    peaks = {}
    for idx in filtered:
        freq_min, freq_max = find_peak_bounds(detrended, signal_freqs, idx, window=31)
        peaks[signal_freqs[idx]] = (freq_min, freq_max)
    return peaks
