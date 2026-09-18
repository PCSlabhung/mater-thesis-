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

A minimum

