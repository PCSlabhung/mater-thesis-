# MATLAB Simulation and Analysis

This branch contains the MATLAB scripts used for air-coupled ultrasonic simulation, synthetic-aperture beamforming, hardware-output visualization, and diagnostic analysis.

## Overview

The current system is based on a 5-Tx / 64-Rx air-ultrasound imaging architecture.

Main parameters:

```matlab
f0 = 40e3;                  % Center frequency [Hz]
c  = 340;                   % Speed of sound [m/s]
fs = 3.125e6 / 19;          % Sampling frequency

no_lines     = 69;          % Azimuth scan lines
no_lines_phi = 69;          % Elevation scan lines

d_th     = 0.025574095887433;
d_th_phi = 0.025574095887433;

range = [300:6:2000] / fs * c / 2;
```

The MATLAB scripts are mainly divided into three parts:

1. Field II simulation and baseband generation
2. MATLAB synthetic-aperture beamforming
3. Hardware-output analysis and diagnostics

---

## File Structure

```text
matlab/
├── STA_64ch_sunflowerV2.m
├── create_bb_dataV2.m
├── hw_output_figure.m
├── beampattern.m
├── SA_analyze.m
├── SA_depth_analyze.m
└── README.md
```

Hardware-output datasets are expected to contain files such as:

```text
dataset/
├── SA_image_sum.txt
├── SA_image_tx01.txt
├── SA_image_tx02.txt
├── SA_image_tx03.txt
├── SA_image_tx04.txt
└── SA_image_tx05.txt
```

---

# 1. create_bb_dataV2.m

## Purpose

`create_bb_dataV2.m` generates the simulated ultrasonic receive data used by the MATLAB beamformer.

The script uses Field II to:

* Construct the Tx aperture
* Construct the 64-channel Rx array
* Define the acoustic parameters
* Simulate point-target echoes
* Generate RF data
* Downsample the RF data
* Perform complex demodulation
* Generate complex baseband data

The resulting data are stored in:

```matlab
d(:,:,tx_index)
```

with the logical format:

```text
d(depth_sample, rx_channel, tx_index)
```

---

## Tx Geometry

Five transmitters are used.

Their positions are defined by:

```matlab
temp = 0.0783 + 15/1000;

right_x = temp;
left_x  = -temp;

top_y = temp;
bot_y = -temp;

vs_x = [0 left_x right_x 0 0];
vs_y = [0 0 0 top_y bot_y];
vs_z = 0;
```

The geometry is approximately:

```text
             Tx4
              |
              |
Tx2 -------- Tx1 -------- Tx3
              |
              |
             Tx5
```

Tx1 is located at the center.

Tx2 and Tx3 are placed along the X direction.

Tx4 and Tx5 are placed along the Y direction.

---

## Rx Geometry

The receiver contains 64 channels arranged using a sunflower / Fermat-spiral distribution.

The aperture radius is:

```matlab
Rap = 0.118 / 2;
```

The radial position is:

```matlab
RR(m) = Rap * sqrt(m/64);
```

The angular position is:

```matlab
phi(m) = 2*pi*(m-1)*(1+sqrt(5))/2;
```

The Cartesian coordinates are:

```matlab
chx(m) = RR(m) * cos(phi(m));
chy(m) = RR(m) * sin(phi(m));
```

A minimum radius is applied to the elements near the array center.

---

## Acoustic Parameters

The current simulation parameters are:

```matlab
f0 = 40e3;
c  = 340;

fs = 3.125e6 / 19;
fs_simu = fs * 4;
```

where:

```text
f0      = ultrasonic center frequency
c       = speed of sound in air
fs      = final sampling frequency
fs_simu = Field II simulation sampling frequency
```

---

## Tx Excitation

The current Tx excitation waveform is generated using PWM.

```matlab
pwm_t = 0:1/fs_simu:50e-6;

pwm1 = double( ...
    mod(pwm_t, 2/80e3) < 1/80e3 ...
);
```

The Tx impulse response is based on the 40 kHz center frequency.

---

## Rx Impulse Response

The Rx impulse response is currently modeled using:

```matlab
impulse_response_RX = ...
    cos(2*pi*28000*(0:1/fs_simu:2.5/f0));
```

---

## Field II Dependency

This script requires Field II.

The current script contains a machine-specific path:

```matlab
path(path, ...
'C:\DSP計畫\DSP計畫\3D\polar\0318\0318\fieldII');
```

Before running on another computer, modify the path so that MATLAB can find Field II.

The script initializes Field II using:

```matlab
field_init
```

---

# 2. STA_64ch_sunflowerV2.m

## Purpose

`STA_64ch_sunflowerV2.m` is the main MATLAB reference beamforming implementation.

It performs:

* 5-Tx synthetic-aperture processing
* 64-channel receive beamforming
* 3-D scan-line reconstruction
* Exact delay calculation
* Hardware-approximate delay calculation
* Delay approximation error analysis
* Coherence Factor processing
* Tx accumulation
* B-mode visualization
* C-mode visualization
* Beampattern analysis
* 3-D imaging visualization

The script calls:

```matlab
create_bb_dataV2
```

to generate the input baseband data.

---

## Imaging Grid

The current scan grid is:

```matlab
no_lines = 69;
no_lines_phi = 69;
```

The angular steps are:

```matlab
d_th = 0.025574095887433;
d_th_phi = 0.025574095887433;
```

Azimuth sampling:

```matlab
sin_th = ...
    -(no_lines-1)/2*d_th ...
    + (jj-1)*d_th;
```

Elevation sampling:

```matlab
sin_phi = ...
    -(no_lines_phi-1)/2*d_th_phi ...
    + (ms-1)*d_th_phi;
```

---

## Imaging Coordinates

For each depth and scan direction, the imaging position is generated as:

```matlab
x_org = range.' * sin_th * cos_phi;
y_org = range.' * sin_phi;
z_org = range.' * cos_th * cos_phi;
```

where:

```matlab
cos_th  = cos(asin(sin_th));
cos_phi = cos(asin(sin_phi));
```

---

# Delay Calculation

Two delay models are implemented.

## Software Delay

The exact propagation distance is calculated using the real Tx and Rx coordinates.

Tx distance:

```matlab
Tx_dist = sqrt( ...
    (x_org - vs_x).^2 + ...
    (z_org - vs_z).^2 + ...
    (y_org - vs_y).^2 ...
);
```

Rx distance:

```matlab
onewaydist = sqrt( ...
    (x_org - chx).^2 + ...
    (z_org - chz).^2 + ...
    (y_org - chy).^2 ...
);
```

The exact delay is:

```matlab
totdelay = ...
    (onewaydist + Tx_dist) / c - toff;
```

---

## Hardware Delay Approximation

The hardware approximation separates the Rx distance into two terms.

The first term is:

```matlab
T_rd = sqrt( ...
    range.^2 + ...
    chx.^2 + ...
    chy.^2 ...
);
```

The direction-dependent correction is:

```matlab
T_sd_rd = -( ...
    chx * cos_phi * sin_th + ...
    chy * sin_phi ...
);
```

The approximate total delay is:

```matlab
totdelay = ...
    (Tx_dist + T_rd + T_sd_rd) / c;
```

The delay model can be selected using:

```matlab
index_mode = 'hardware';
```

or:

```matlab
index_mode = 'software';
```

---

# Delay Approximation Error Analysis

The script compares:

```text
Exact Rx distance
vs.
Hardware approximate Rx distance
```

The following errors are calculated:

```text
Distance error [mm]
Delay error [ns]
Sample error [samples]
Phase error [degree]
```

For each error type, the script reports:

```text
Mean absolute error
RMSE
Maximum absolute error
```

The script also records the location of the maximum sample error, including:

```text
Tx index
Theta index
Phi index
Theta angle
Phi angle
Depth index
Depth
Rx channel
```

This part is useful for evaluating the effect of the hardware delay approximation.

---

# Beamforming

The calculated delay is converted to a sample index:

```matlab
delayindx = totdelay * fs + 1;
```

The corresponding baseband sample is selected from each Rx channel.

Phase compensation is then applied:

```matlab
exp(j*2*pi*f0*totdelay)
```

The channel data are:

```matlab
chandata = ...
    interp_bbdata ...
    .* exp(j*2*pi*f0*totdelay);
```

The DAS beamforming result is:

```matlab
subbf(:,jj,ms) = ...
    sum(chandata, 2);
```

---

# Coherence Factor

The script calculates a Coherence Factor using:

```matlab
CF = ...
    abs(subbf).^2 ...
    ./ ...
    (M * sum(abs(chandata.^2),2)) ...
    + 1e-8;
```

The beamforming output is then weighted using:

```matlab
CF.^2
```

---

# PAF

The current PAF parameter is:

```matlab
alph = 0;
```

The implemented expression is:

```matlab
subbf_paf = ...
    (1-alph)*subbf ...
    + alph*(subbf.^2);
```

Since:

```matlab
alph = 0;
```

the nonlinear PAF term is currently disabled.

---

# Synthetic Aperture Accumulation

The output from the five Tx positions is accumulated using:

```matlab
all_subbf = ...
    all_subbf ...
    + subbf_paf .* CF.^2;
```

The final `all_subbf` contains the synthetic-aperture reconstruction.

---

# MATLAB Visualization

After beamforming, the script performs log compression:

```matlab
bfn_data = ...
    abs(all_subbf / max(all_subbf(:)));

qqq = 20*log10(bfn_data);
```

The default dynamic range is:

```matlab
DR = 40;
```

---

## XZ B-mode

The XZ image is generated using:

```matlab
X = range.' * sin_th;
Z = range.' * cos_th;
```

and displayed using:

```matlab
pcolor(X, Z, QQ)
```

---

## XY C-mode

A depth region around the detected target is selected.

The script then performs a maximum-intensity projection across several depth samples.

The result is interpolated to a higher resolution for visualization.

---

## Beampattern

The script extracts the azimuth beampattern near the target.

The beam is interpolated using:

```matlab
interp1(..., 'spline')
```

The `-6 dB` width is calculated from the interpolated curve.

---

## 3-D Image

The script performs 3-D interpolation from:

```text
Range
Azimuth
Elevation
```

to Cartesian coordinates.

The result is visualized using MATLAB `slice`.

---

# 3. hw_output_figure.m

## Purpose

`hw_output_figure.m` is used to visualize one reconstructed hardware-output file.

The current input file is manually selected using:

```matlab
filename = ...
    '30_30_100_16bit/SA_image_sum.txt';
```

Change this path before analyzing another dataset.

---

## Expected Input Format

The hardware TXT file must contain:

```text
rows    = no_lines × no_lines_phi
columns = Nr
```

For the current configuration:

```text
rows = 69 × 69 = 4761
```

The data are converted using:

```matlab
bf = reshape( ...
    loaded_data.', ...
    [Nr, no_lines, no_lines_phi] ...
);
```

Therefore the internal data format is:

```text
bf(depth, azimuth, elevation)
```

---

## Automatic Target Detection

The script automatically searches for the strongest voxel:

```matlab
[max_val, max_idx_3d] = max(qqq(:));
```

The target indices are:

```matlab
[idx_r, idx_az, idx_el]
```

The corresponding physical location is estimated as:

```text
Range
Azimuth
Elevation
X
Y
Z
```

---

## XZ B-mode

The script automatically selects the elevation plane containing the strongest target.

The result is displayed as:

```text
Azimuth X
vs.
Depth Z
```

---

## XY C-mode

A small depth region around the detected target is used.

The script performs Maximum Intensity Projection across the depth region.

The C-mode image therefore represents:

```text
Azimuth X
vs.
Elevation Y
```

near the target depth.

---

## Beampattern

The script uses a small depth slab around the target:

```matlab
r_margin_bp = 5;
```

and performs:

```matlab
bp_raw = max(slab, [], 1);
```

This avoids accidentally selecting an empty depth sample.

The extracted beampattern is interpolated using spline interpolation.

---

## -6 dB Beamwidth

The threshold is defined as:

```matlab
threshold_dB = peak_dB - 6;
```

The script searches left and right from the main peak to estimate the main-lobe width.

The width is reported in meters.

---

## 3-D Volume Slice

The script generates a local 3-D volume around the detected target.

The current zoom range is:

```matlab
zoom = 0.2;
```

which corresponds to approximately:

```text
±20 cm
```

around the detected target.

---

# 4. beampattern.m

## Purpose

`beampattern.m` compares the reconstructed beampatterns from:

```text
Tx01
Tx02
Tx03
Tx04
Tx05
SA Sum
```

The current file paths are defined manually:

```matlab
filename1 = ...
    '0_0_100_12bit/SA_image_tx01.txt';

filename2 = ...
    '0_0_100_12bit/SA_image_tx02.txt';

filename3 = ...
    '0_0_100_12bit/SA_image_tx03.txt';

filename4 = ...
    '0_0_100_12bit/SA_image_tx04.txt';

filename5 = ...
    '0_0_100_12bit/SA_image_tx05.txt';

filename_sum = ...
    '0_0_100_12bit/SA_image_sum.txt';
```

Change these paths before analyzing another dataset.

---

## Processing

For every Tx file and the SA Sum file, the script:

1. Loads the TXT data
2. Checks the data dimensions
3. Reshapes the data into a 3-D image
4. Normalizes the image
5. Converts the result to dB
6. Automatically locates the strongest target
7. Selects a small depth slab around the target
8. Performs a depth MIP
9. Extracts the azimuth beampattern
10. Performs spline interpolation
11. Calculates the `-6 dB` beamwidth

---

## Fair Beamwidth Comparison

Each beampattern uses its own peak as the reference.

For example:

```text
Tx01 threshold = Tx01 peak - 6 dB
Tx02 threshold = Tx02 peak - 6 dB
...
Sum threshold = Sum peak - 6 dB
```

This avoids comparing different curves using an incorrect common absolute threshold.

---

## Output Figures

The script generates:

### Absolute Beampattern Comparison

Displays:

```text
Tx01
Tx02
Tx03
Tx04
Tx05
Sum
```

on the same plot.

---

### Normalized Beampattern Comparison

Each beampattern is shifted so that its own maximum is:

```text
0 dB
```

A common:

```text
-6 dB
```

reference line is then used for shape comparison.

---

## Reported Results

The command window displays:

```text
Peak amplitude
-6 dB beamwidth
```

for:

```text
Tx01
Tx02
Tx03
Tx04
Tx05
Sum
```

---

# 5. SA_analyze.m

## Purpose

`SA_analyze.m` is a batch diagnostic script used to compare:

```text
Individual Tx beamwidth
vs.
Synthetic Aperture Sum beamwidth
```

It is mainly intended for debugging cases where the SA image is not narrower than the individual Tx images.

---

## Input Folder

The current folder is:

```matlab
folder_path = ...
    '30_30_150_float/';
```

The script automatically searches for:

```matlab
SA_image_tx*.txt
```

---

## Processing

The script first processes:

```text
SA_image_sum.txt
```

and detects:

```text
SA Sum peak depth
SA Sum -6 dB beamwidth
```

It then processes every:

```text
SA_image_txXX.txt
```

file in the same folder.

---

## Depth Filtering

A Tx result is only included if its detected peak depth is close to the SA Sum target.

The current criterion is:

```matlab
abs(d - sum_peak_depth) < 0.2
```

which corresponds to:

```text
±0.2 m
```

This prevents obvious noise peaks from being included in the comparison.

---

## Results

The script calculates:

```text
SA Sum beamwidth
Average single-Tx beamwidth
Minimum single-Tx beamwidth
```

It also plots:

```text
Tx Index
vs.
-6 dB Beamwidth
```

with the SA Sum beamwidth shown as a reference line.

---

# 6. SA_depth_analyze.m

## Purpose

`SA_depth_analyze.m` checks whether different Tx reconstructions are aligned in the depth direction.

This is important because depth misalignment between Tx images can cause:

```text
Destructive summation
Reduced synthetic-aperture gain
Beam broadening
Defocusing
Target position shift
```

---

## Input Folder

The current folder is:

```matlab
folder_path = '0_0_150/';
```

---

## Target Locking

The script first loads:

```text
SA_image_sum.txt
```

and finds the strongest target.

The detected:

```text
Azimuth index
Elevation index
```

are then fixed for every Tx.

This means all A-scans are compared along exactly the same imaging direction.

---

## Depth Profile Comparison

For every Tx:

```text
Tx01
Tx02
Tx03
Tx04
Tx05
```

the script extracts:

```matlab
profile = ...
    data(:, idx_az, idx_el);
```

The profile is normalized and converted into dB.

All five depth profiles are displayed on the same figure.

---

## Tx01 vs Tx05 Diagnosis

The second diagnostic figure specifically compares:

```text
Tx01
vs.
Tx05
```

The peak depth of both Tx images is detected.

The script calculates:

```matlab
diff_m = depth01 - depth05;
```

and reports the difference in:

```text
meters
centimeters
```

This comparison can help identify Tx-dependent delay or indexing errors.

---

# Hardware Output Data Format

The MATLAB hardware-analysis scripts expect the hardware reconstruction output to have the following format:

```text
Number of rows:
no_lines × no_lines_phi

Number of columns:
Nr
```

With the current parameters:

```text
69 × 69 = 4761 rows
```

The conversion is:

```matlab
volume = reshape( ...
    raw.', ...
    [Nr, no_lines, no_lines_phi] ...
);
```

The internal representation is:

```text
volume(depth, azimuth, elevation)
```

---

# Recommended Workflow

## 1. MATLAB Reference Simulation

Run:

```matlab
STA_64ch_sunflowerV2
```

This automatically runs:

```matlab
create_bb_dataV2
```

Use this flow to verify:

```text
Tx/Rx geometry
Field II simulation
Delay equations
Hardware delay approximation
DAS beamforming
Coherence Factor
Synthetic-aperture accumulation
B-mode
C-mode
Beampattern
```

---

## 2. Analyze Hardware SA Output

Modify:

```matlab
filename
```

inside:

```text
hw_output_figure.m
```

Then run:

```matlab
hw_output_figure
```

This gives a general visualization of:

```text
SA_image_sum.txt
```

including:

```text
Target location
B-mode
C-mode
Beampattern
-6 dB beamwidth
3-D slice
```

---

## 3. Compare Tx Images and SA Sum

Modify the file paths inside:

```text
beampattern.m
```

Then run:

```matlab
beampattern
```

Use this script when comparing the spatial resolution of:

```text
Single Tx reconstruction
vs.
Synthetic Aperture reconstruction
```

---

## 4. Diagnose SA Beamwidth

If the synthetic-aperture beamwidth appears incorrect, run:

```matlab
SA_analyze
```

This shows whether:

```text
SA Sum beamwidth
```

is smaller or larger than the individual Tx beamwidths.

---

## 5. Diagnose Depth Alignment

If the SA Sum appears defocused or wider than expected, run:

```matlab
SA_depth_analyze
```

This checks whether the Tx images are aligned in depth before accumulation.

---

# Parameters That Must Match the Hardware

The following parameters are repeated in multiple MATLAB scripts and must match the FPGA/HLS configuration.

```matlab
c = 340;

fs = 3.125e6 / 19;

no_lines = 69;

no_lines_phi = 69;

d_th = 0.025574095887433;

d_th_phi = 0.025574095887433;

range = ...
    [300:6:2000] / fs * c / 2;
```

Also verify:

```text
Tx coordinates
Rx coordinates
Rx channel ordering
Number of Rx channels
Number of Tx
Depth sample ordering
Azimuth ordering
Elevation ordering
Fixed-point scaling
Output TXT ordering
```

Any mismatch between MATLAB and hardware may cause:

```text
Incorrect target position
Incorrect range
Incorrect beampattern
Incorrect beamwidth
Tx depth misalignment
Synthetic-aperture defocusing
```

---

# Current Caveats

## Hard-Coded Dataset Paths

Several scripts contain experiment-specific paths.

Examples:

```matlab
filename = ...
    '30_30_100_16bit/SA_image_sum.txt';
```

```matlab
folder_path = ...
    '30_30_150_float/';
```

```matlab
folder_path = ...
    '0_0_150/';
```

These paths must be changed manually when using another dataset.

---

## Field II Path

`create_bb_dataV2.m` contains a local Windows Field II path.

This should eventually be replaced by either:

```text
Relative path
```

or:

```text
MATLAB setup script
```

---

## Shared Workspace

`STA_64ch_sunflowerV2.m` currently calls:

```matlab
create_bb_dataV2
```

as a script.

Therefore both scripts share the same MATLAB workspace.

Variables generated inside `create_bb_dataV2.m`, such as:

```matlab
d
```

are directly used by `STA_64ch_sunflowerV2.m`.

Because of this, moving the files into separate folders may require adding the folders to the MATLAB path.

---

## Repeated Parameters

Parameters such as:

```text
fs
c
no_lines
no_lines_phi
d_th
range
```

are currently repeated in several scripts.

If one script is modified while another is not, MATLAB and hardware analysis may become inconsistent.

A future improvement would be to move all shared parameters into one configuration file.

For example:

```text
config.m
```

---

# Suggested Future Structure

A cleaner folder structure could be:

```text
matlab/
├── config/
│   └── imaging_config.m
│
├── simulation/
│   ├── create_bb_dataV2.m
│   └── STA_64ch_sunflowerV2.m
│
├── analysis/
│   ├── hw_output_figure.m
│   ├── beampattern.m
│   ├── SA_analyze.m
│   └── SA_depth_analyze.m
│
└── README.md
```

However, the current flat structure is easier to use because:

```matlab
STA_64ch_sunflowerV2
```

directly calls:

```matlab
create_bb_dataV2
```

as a script.

---

# Script Summary

| Script                   | Purpose                                             |
| ------------------------ | --------------------------------------------------- |
| `create_bb_dataV2.m`     | Field II simulation and complex baseband generation |
| `STA_64ch_sunflowerV2.m` | MATLAB STA + CF beamforming reference               |
| `hw_output_figure.m`     | Hardware SA output visualization                    |
| `beampattern.m`          | Tx01-Tx05 and SA Sum beampattern comparison         |
| `SA_analyze.m`           | Batch Tx vs SA beamwidth diagnostic                 |
| `SA_depth_analyze.m`     | Tx depth-alignment diagnostic                       |

---

# Typical Analysis Flow

```text
Field II
   │
   ▼
create_bb_dataV2.m
   │
   │ complex baseband data
   ▼
STA_64ch_sunflowerV2.m
   │
   │ MATLAB reference image
   ▼
Software / Hardware Comparison
```

For hardware output:

```text
FPGA / HLS Output
   │
   ├── SA_image_tx01.txt
   ├── SA_image_tx02.txt
   ├── SA_image_tx03.txt
   ├── SA_image_tx04.txt
   ├── SA_image_tx05.txt
   └── SA_image_sum.txt
            │
            ▼
      MATLAB Analysis
            │
            ├── hw_output_figure.m
            ├── beampattern.m
            ├── SA_analyze.m
            └── SA_depth_analyze.m
```

---

# Summary

This MATLAB branch provides the software reference and analysis environment for the 5-Tx / 64-Rx air-coupled ultrasonic synthetic-aperture imaging system.

It supports:

```text
Field II simulation
64-channel sunflower Rx array
5-Tx synthetic aperture
3-D DAS beamforming
Coherence Factor processing
Hardware delay approximation verification
B-mode imaging
C-mode imaging
Beampattern analysis
-6 dB beamwidth measurement
Tx-to-Tx depth alignment analysis
Hardware-output visualization
```

The scripts are primarily intended for verifying the FPGA/HLS beamforming implementation and diagnosing differences between MATLAB and hardware reconstruction results.

