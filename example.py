import time
import numpy as np
import matplotlib.pyplot as plt
from collections import deque

from signal_generator.signal_generator import RandomSignalGenerator
from peak_detector.peak_detector import (
    detect_peaks,
    _calculate_wavelet_coeffs,
    _find_local_maxima,
)

STEP_FREQUENCY = 124
EQ_RANGE = 10000
MAX_ITERATIONS = 1000


def plot_signal(ax, x, y, detected_peaks):
    ax.cla()
    ax.plot(x, y, label="Signal")
    for freq in detected_peaks:
        ax.axvline(x=freq, color="r", linestyle="--", label="Detected peak")
    ax.set_title("Real-time Peak Detection")
    ax.set_xlabel("Frequency [Hz]")
    ax.set_ylabel("Amplitude")
    ax.legend(loc="upper right")
    ax.grid(True)


def show_debug_plot(x, y, detected_peaks):
    coeffs = _calculate_wavelet_coeffs(np.array(y), STEP_FREQUENCY)
    maxima = _find_local_maxima(coeffs)

    plt.ioff()
    fig, axs = plt.subplots(2, figsize=(10, 6))
    axs[0].plot(x, y, label="Signal")
    for freq in detected_peaks:
        axs[0].axvline(x=freq, color="r", linestyle="--")
    axs[0].set_title("Failed Detection - Signal")
    axs[0].legend()

    axs[1].imshow(maxima, aspect="auto")
    axs[1].set_title("Wavelet Coefficients (Local Maxima)")
    plt.tight_layout()
    plt.show()


def main():
    generator = RandomSignalGenerator()

    plt.ion()
    fig, ax = plt.subplots(figsize=(10, 4))

    try:
        for _ in range(MAX_ITERATIONS):
            x = deque(
                np.arange(
                    generator.FREQUENCY_RANGE[0],
                    generator.FREQUENCY_RANGE[1],
                    STEP_FREQUENCY,
                )
            )
            y = deque(generator.generate_random_signal(val) for val in x)

            start = time.time()
            detected = detect_peaks(x, y)
            duration = time.time() - start
            print(f"Detect time: {duration:.4f}s")

            if len(generator.peaks) != len(detected):
                raise ValueError("Number of peaks mismatch")

            for freq in detected:
                min_diff = min(
                    abs(freq - peak["frequency"]) for peak in generator.peaks
                )
                if min_diff > EQ_RANGE:
                    raise ValueError(f"Frequency mismatch: diff={min_diff}")

            plot_signal(ax, x, y, detected)
            plt.pause(0.01)
            generator.generate_new_number_of_peaks()
            generator.generate_new_signal_params()

    except Exception as e:
        print(f"[ERROR] {e}")
        show_debug_plot(x, y, list(detected))


if __name__ == "__main__":
    main()
