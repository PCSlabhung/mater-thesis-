# Test Data

This directory contains the input datasets used to verify the synthetic aperture (SA) ultrasound beamforming implementation.

Each dataset represents a point target located at a specific angular position and range.

## Directory Structure

```text
TEST_DATA/
├── -30_-30_50/
├── -30_-30_100/
├── -30_-30_150/
│
├── -30_30_50/
├── -30_30_100/
├── -30_30_150/
│
├── 0_0_50/
├── 0_0_100/
├── 0_0_150/
│
├── 30_-30_50/
├── 30_-30_100/
├── 30_-30_150/
│
├── 30_30_50/
├── 30_30_100/
└── 30_30_150/
```

## Dataset Naming Convention

Each dataset directory follows the naming convention:

```text
<theta>_<phi>_<range>
```

where:

| Field | Description |
|---|---|
| `theta` | Target angle in the first scanning dimension |
| `phi` | Target angle in the second scanning dimension |
| `range` | Target distance in centimeters |

For example:

```text
0_0_100
```

represents a target located at:

```text
theta = 0°
phi   = 0°
range = 100 cm
```

Similarly:

```text
30_-30_150
```

represents:

```text
theta = 30°
phi   = -30°
range = 150 cm
```

The datasets currently cover three target ranges:

```text
50 cm
100 cm
150 cm
```

and multiple angular positions:

```text
(-30°, -30°)
(-30°,  30°)
(  0°,   0°)
( 30°, -30°)
( 30°,  30°)
```

---

## Files in Each Dataset

Each dataset contains the complex baseband receive signals for five Tx positions.

For example:

```text
0_0_100/
├── bb_tx01_real.txt
├── bb_tx01_imag.txt
├── bb_tx02_real.txt
├── bb_tx02_imag.txt
├── bb_tx03_real.txt
├── bb_tx03_imag.txt
├── bb_tx04_real.txt
├── bb_tx04_imag.txt
├── bb_tx05_real.txt
└── bb_tx05_imag.txt
```

Each transmitter contains two files:

```text
bb_txXX_real.txt
bb_txXX_imag.txt
```

where `XX` is the Tx index from `01` to `05`.

For example:

```text
bb_tx01_real.txt
bb_tx01_imag.txt
```

contain the real and imaginary components of the complex baseband receive data for Tx 1.

---

## Transmitter Index

The five Tx datasets correspond to the five transmitter positions used in the synthetic aperture imaging system:

```text
           Tx4
            |
            |
Tx3 ------ Tx1 ------ Tx2
            |
            |
           Tx5
```

The Tx numbering is:

| Tx | Position |
|---|---|
| Tx1 | Center |
| Tx2 | +x |
| Tx3 | -x |
| Tx4 | +y |
| Tx5 | -y |

The five Tx measurements are processed independently by the beamformer and subsequently accumulated to generate the synthetic aperture image.

---

## Receive Data

Each Tx file contains receive data from the 64-channel ultrasonic receiver array.

The beamforming implementation uses:

```text
64 receive channels
2000 samples per channel
```

The real and imaginary components are stored separately.

Conceptually, the input data correspond to:

```text
Tx
 │
 └── 64 Rx channels
       │
       ├── Channel 1  ── 2000 samples
       ├── Channel 2  ── 2000 samples
       ├── ...
       └── Channel 64 ── 2000 samples
```

These data are loaded by the Vitis HLS testbench and used as the input to the beamforming architecture.

---

## Synthetic Aperture Processing

For each dataset, the five Tx measurements are processed sequentially:

```text
Tx1 ──► Beamformer ──► Partial Image 1
Tx2 ──► Beamformer ──► Partial Image 2
Tx3 ──► Beamformer ──► Partial Image 3
Tx4 ──► Beamformer ──► Partial Image 4
Tx5 ──► Beamformer ──► Partial Image 5
                         │
                         ▼
                Synthetic Aperture
                   Accumulation
                         │
                         ▼
                    Final Image
```

The resulting partial images are accumulated to form the final synthetic aperture image.

---

## Example

To test a target located at:

```text
theta = 30°
phi   = 30°
range = 100 cm
```

use the dataset:

```text
30_30_100/
```

The corresponding input files are:

```text
30_30_100/
├── bb_tx01_real.txt
├── bb_tx01_imag.txt
├── bb_tx02_real.txt
├── bb_tx02_imag.txt
├── bb_tx03_real.txt
├── bb_tx03_imag.txt
├── bb_tx04_real.txt
├── bb_tx04_imag.txt
├── bb_tx05_real.txt
└── bb_tx05_imag.txt
```

---

## Testbench Usage

The dataset used by the C++ testbench can be selected in `SA_v3_tb.cpp`.

For example:

```cpp
vector<string> datasets = {
    "0_0_100"
};
```

Multiple datasets can also be tested sequentially:

```cpp
vector<string> datasets = {
    "-30_-30_50",
    "-30_-30_100",
    "-30_-30_150",

    "-30_30_50",
    "-30_30_100",
    "-30_30_150",

    "0_0_50",
    "0_0_100",
    "0_0_150",

    "30_-30_50",
    "30_-30_100",
    "30_-30_150",

    "30_30_50",
    "30_30_100",
    "30_30_150"
};
```

---

## Notes

These datasets are intended for functional verification of the FPGA/HLS synthetic aperture ultrasound beamformer.

They are used to evaluate beamforming results for different target directions and ranges and to verify the depth-sparse synthetic aperture accumulation architecture.
