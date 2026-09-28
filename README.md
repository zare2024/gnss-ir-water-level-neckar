# GNSS-IR water level: Neckar river, Ober-Esslingen

Estimating the water level of the Neckar from GNSS signals reflected off the river surface (GNSS interferometric reflectometry). A receiver on the riverbank records the signal-to-noise ratio (SNR) of GPS, Galileo and GLONASS satellites. The interference between direct and reflected signals makes the SNR oscillate, and the frequency of that oscillation gives the height of the antenna above the water.

**Result: reflector height Hr = 3.04 ± 0.18 m** (8 consistent satellite arcs).

![Reflector height of the kept arcs](figures/arc_heights.png)

## Method

For a horizontal reflecting surface, the reflected part of the SNR follows

```
SNR_ref = A · sin(2π f · sin(ε) + φ),    f = 2 Hr / λ
```

where ε is the satellite elevation angle, Hr the reflector height and λ the carrier wavelength (L1, 0.1903 m).

`gnss_ir_water_level.m` processes every satellite in the dataset:

1. **Azimuth sector, 270° to 330°.** Keeps reflections from the river and drops those from land.
2. **Arc split.** Sorts each satellite's samples by time and splits the pass at its highest elevation into a rising and a setting arc, which are processed separately.
3. **Direct-signal removal.** Fits a 2nd-order polynomial to SNR as a function of sin(ε) and subtracts it.
4. **Spectral analysis.** Computes a Lomb–Scargle periodogram (samples are unevenly spaced in sin(ε)) and converts the dominant frequency to Hr = f · λ / 2. Arcs with fewer than 100 samples are skipped.
5. **Robust combination.** Keeps arcs within ±1.5 m of the median and reports their mean and standard deviation, overall and per constellation.

All parameters sit in the settings block at the top of the script.

## Results

| | |
|---|---|
| Median of all arcs | 3.10 m |
| Arcs kept | 8 |
| Outlier arcs removed | e.g. 124 m, 203 m, 500 m |
| **Reflector height** | **3.04 ± 0.18 m** |

Estimates from GPS, Galileo and GLONASS were mostly consistent; the differences were confined to the outlier arcs. The script prints a per-constellation mean to check this.

## Known limitations and next steps

- **Frequency grid.** The periodogram searches Hr from 0.1 to 500 m with 2,000 frequencies, a step of about 0.25 m in Hr. That is coarse compared with the 0.18 m spread. A narrower search range (e.g. 0.5 to 10 m) would give millimetre-level steps and would also avoid the spurious peaks at hundreds of metres.
- **GLONASS wavelength.** GLONASS transmits on slightly different L1 frequencies per satellite (λ ≈ 0.187 m), but the script uses the GPS L1 wavelength for all satellites. This scales GLONASS heights by about 1 to 2 %.
- **Elevation mask.** Limiting the arcs to low-to-mid elevations (e.g. 5° to 25°) should strengthen the reflection signal.
- **Time series.** Estimating Hr in a moving window (e.g. 30 to 60 min) would turn the single estimate into a water level time series.

## Data

The SNR observations are not included. They were recorded by the Institute of Geodesy, University of Stuttgart, on 31 October 2022 and provided for the course.

The script expects `gnss_ir.mat` next to it: a MATLAB struct `gnss_ir` with one table per satellite (10 Galileo `E..`, 15 GPS `G..`, 12 GLONASS `R..`), each with the variables `time`, `azimut`, `elevation` and `snr1` (L1 SNR).

## How to run

Requires MATLAB R2016b or newer with the Signal Processing Toolbox (`plomb`).

```matlab
gnss_ir_water_level
```

The script prints one line per satellite arc, the combined result and a per-constellation summary, and saves `figures/arc_heights.png`.

## Credits

Author: Shaghayegh Zare. Developed from Lab 3 of the course *Satellite Geodesy Observation Techniques* (M.Sc. Geomatics Engineering, University of Stuttgart, winter semester 2025/26), completed in a team with Haibo and Shiwei; the code was rewritten and restructured for this repository.

References: Larson et al. (2009), IEEE JSTARS 3(1), 91–99; Lomb (1976), Astrophysics and Space Science 39(2), 447–462; Scargle (1982), The Astrophysical Journal 263, 835–853.
