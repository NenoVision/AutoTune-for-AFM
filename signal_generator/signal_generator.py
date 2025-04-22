import numpy as np
import random
import scipy


class RandomSignalGenerator:
    FREQUENCY_RANGE = (2000, 500000)
    PEAK_COUNT_RANGE = (1, 3)
    NOISE_LEVEL = 0.002

    PEAK_PARAMS = {
        "width": (10_000, 30_000),
        "height": (0.04, 1.1),
        "shannonCoeff": (1, 10),
        "shannonOffset": (np.pi / 2, np.pi / 10),
        "poissonCoeff": (5, 30),
        "mexicanHatCoeff": (0.01, 0.25),
    }

    SECOND_PEAK_FREQ = (24_000, 82_000)
    THIRD_PEAK_FREQ = (135_000, 275_800)

    def __init__(self, number_of_peaks: int = 0) -> None:
        try:
            self.number_of_peaks = number_of_peaks
            self._validate_peak_count()
        except ValueError:
            self.generate_new_number_of_peaks()

        self.generate_new_signal_params()

    def generate_new_number_of_peaks(self) -> None:
        self.number_of_peaks = self._randint(*self.PEAK_COUNT_RANGE)

    def _validate_peak_count(self):
        min_peaks, max_peaks = self.PEAK_COUNT_RANGE
        if not (min_peaks <= self.number_of_peaks <= max_peaks):
            raise ValueError("Number of peaks is out of the range.")

    def generate_new_signal_params(self):
        self.bezier_points = self._generate_bezier_points()
        self.peaks = self._generate_peaks()
        self.noise = self.NOISE_LEVEL
        self.noise_peaks = self._generate_noise_peaks()

    @staticmethod
    def _rand(a: float, b: float, round_val: bool = False) -> float:
        val = random.uniform(a, b)
        return round(val) if round_val else val

    @staticmethod
    def _randint(a: int, b: int) -> int:
        return random.randint(a, b)

    def _generate_noise_peaks(self):
        return [self._create_random_peak(height_range=(0.003, 0.01)) for _ in range(4)]

    def _create_random_peak(self, height_range=None):
        peak = {k: self._rand(*v) for k, v in self.PEAK_PARAMS.items()}
        peak["frequency"] = self._rand(*self.FREQUENCY_RANGE)
        if height_range:
            peak["height"] = self._rand(*height_range)
        return peak

    def _generate_bezier_points(self):
        start_val = 0
        end_val = self._rand(0, 0.2)
        peak_freq = self._rand(300_000, 450_000)
        peak_val = self._rand(0, 0.2)

        lower_freq = peak_freq - self._rand(10_000, 80_000)
        upper_freq = peak_freq + self._rand(10_000, 40_000)

        lower_val = peak_val + self._rand(0, 0.3)
        upper_val = peak_val + self._rand(0, 0.3)

        return np.array(
            [
                [self.FREQUENCY_RANGE[0], start_val],
                [lower_freq, lower_val],
                [peak_freq, peak_val],
                [upper_freq, upper_val],
                [self.FREQUENCY_RANGE[1], end_val],
            ]
        )

    def _generate_peaks(self):
        peaks = []
        for i in range(self.number_of_peaks):
            peak = {k: self._rand(*v) for k, v in self.PEAK_PARAMS.items()}
            if i == 0:
                peak["frequency"] = self.bezier_points[2][0]
            elif i == 1:
                peak["frequency"] = self._rand(*self.SECOND_PEAK_FREQ)
            elif i == 2:
                peak["frequency"] = self._rand(*self.THIRD_PEAK_FREQ)
            peaks.append(peak)
        return peaks

    def _bezier_value(self, freq: int) -> float:
        t = (freq - self.FREQUENCY_RANGE[0]) / (
            self.FREQUENCY_RANGE[1] - self.FREQUENCY_RANGE[0]
        )
        return self._bezier_curve(t)[1]

    def _bezier_curve(self, t: float) -> float:
        n = len(self.bezier_points) - 1
        return sum(
            scipy.special.comb(n, i) * (1 - t) ** (n - i) * t**i * np.array(pt)
            for i, pt in enumerate(self.bezier_points)
        )

    def _generate_wavelet_sum(self, freq: int, peak: dict) -> float:
        min_f = peak["frequency"] - peak["width"] / 2
        max_f = peak["frequency"] + peak["width"] / 2
        t_scaled = 2 * (freq - min_f) / (max_f - min_f) - 1

        shannon = np.sinc(t_scaled * peak["shannonCoeff"]) * np.cos(
            peak["shannonCoeff"] * np.pi * t_scaled + peak["shannonOffset"]
        )
        poisson = np.exp(-np.abs(t_scaled * peak["poissonCoeff"]))
        mexican = (1 - t_scaled**2 / peak["mexicanHatCoeff"] ** 2) * np.exp(
            -(t_scaled**2) / (2 * peak["mexicanHatCoeff"] ** 2)
        )

        return (shannon + poisson + mexican) / 3

    def generate_random_signal(self, frequency: int) -> float:
        if not (self.FREQUENCY_RANGE[0] <= frequency <= self.FREQUENCY_RANGE[1]):
            raise ValueError("Frequency is out of the range.")

        value = self._bezier_value(frequency)
        for peak in self.peaks + self.noise_peaks:
            value += peak["height"] * self._generate_wavelet_sum(frequency, peak)

        value += self._rand(-self.noise, self.noise)
        return value
