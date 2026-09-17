# Synthetic Aperture Ultrasound Beamformer – Vitis HLS Implementation

This directory contains the C++/Vitis HLS implementation of the synthetic aperture (SA) ultrasound beamforming architecture developed for this thesis.

The design targets real-time air-coupled ultrasonic 3-D imaging and includes delay generation, phase compensation, coherence-factor (CF) beamforming, and depth-sparse synthetic aperture accumulation.

## File Structure

```text
CODE/
├── SA_v3.cpp
├── SA_v3.h
├── SA_v3_tb.cpp
├── sin_cos_value.h
└── README.md
```

## Files

### `SA_v3.cpp`

Main Vitis HLS implementation.

The processing flow is:

```text
Input Real / Imaginary Data
        │
        ▼
Local Channel Cache
        │
        ▼
Beamforming Sequence Generator
        │
        ▼
Delay Approximation
        │
        ▼
Phase Compensation
        │
        ▼
DAS + Coherence Factor Beamformer
        │
        ▼
Depth-Sparse Synthetic Aperture Accumulation
        │
        ▼
Output Image
```

The top-level HLS function is:

```cpp
void top_model(
    bool reset,
    ap_uint<3> TX_idx,
    bf_data_t_int rx_data_real[Data_dim * channel_num],
    bf_data_t_int rx_data_imag[Data_dim * channel_num],
    hls::stream<bf_data_t_64> &out_data,
    hls::stream<bf_data_t_64> &SA_partial_image
);
```

The implementation uses HLS optimization techniques including:

- Dataflow execution
- Loop pipelining
- Loop unrolling
- Array partitioning
- Fixed-point arithmetic
- AXI interfaces
- HLS streams

---

### `SA_v3.h`

Contains the system parameters, fixed-point data types, transducer coordinates, and lookup tables used by the beamformer.

Main system parameters:

| Parameter | Value | Description |
|---|---:|---|
| `channel_num` | 64 | Number of receive channels |
| `NR` | 284 | Number of depth samples |
| `no_lines` | 69 | Number of scan lines |
| `no_lines_phi` | 69 | Number of scan lines in the second angular dimension |
| `Data_dim` | 2000 | Input samples per receive channel |
| `K` | 40 | Number of retained sparse partial sums per Tx |
| Number of Tx | 5 | Number of synthetic aperture transmitters |
| Speed of sound | 340 m/s | Speed of sound in air |
| LUT size | 1024 | Sine/cosine lookup table size |

The transmitter configuration consists of five transmit positions:

```text
        Tx4
         |
Tx3 ---- Tx1 ---- Tx2
         |
        Tx5
```

The receiver array contains 64 elements. Their coordinates are stored in:

```cpp
chx[channel_num]
chy[channel_num]
```

Fixed-point arithmetic is used throughout the design to reduce FPGA hardware resource usage.

---

### `SA_v3_tb.cpp`

C++ testbench used for functional verification of the HLS implementation.

The testbench performs the following operations:

1. Loads the real and imaginary input data.
2. Processes the five transmitters sequentially.
3. Calls the `top_model()` HLS function.
4. Reads the beamformed output stream.
5. Reads the sparse synthetic aperture partial-image stream.
6. Generates individual Tx beamformed images.
7. Accumulates the five Tx results.
8. Checks the reconstructed synthetic aperture image.
9. Generates Top-K and debugging output files.

## Input Data

The input files use the following naming convention:

```text
bb_tx01_real.txt
bb_tx01_imag.txt

bb_tx02_real.txt
bb_tx02_imag.txt

bb_tx03_real.txt
bb_tx03_imag.txt

bb_tx04_real.txt
bb_tx04_imag.txt

bb_tx05_real.txt
bb_tx05_imag.txt
```

Each Tx therefore contains one real-data file and one imaginary-data file.

A recommended repository structure is:

```text
mater-thesis/
├── CODE/
│   ├── SA_v3.cpp
│   ├── SA_v3.h
│   ├── SA_v3_tb.cpp
│   ├── sin_cos_value.h
│   └── README.md
│
└── TEST_DATA/
    └── 0_0_100/
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

Dataset names can be used to indicate the target angular position and distance.

For example:

```text
0_0_100
```

represents one of the test configurations used by the testbench.

---

## Beamforming Architecture

### 1. Beamforming Sequence Generation

The beamformer scans through the image space using three nested loops:

```text
phi
 └── theta
      └── depth
```

For each scan direction and depth sample, the corresponding trigonometric values and range information are generated.

The processing order allows all depth samples belonging to the same scan direction to be processed sequentially, which is useful for the sparse depth-memory architecture.

---

## 2. Delay Approximation

The propagation delay consists of:

```text
Tx → Image Point → Rx
```

The Tx distance is calculated based on the current transmitter position.

For the receiver path, a hardware-efficient approximation is used to reduce the computational cost compared with directly calculating a square root for every receive channel.

The resulting propagation distance is converted to a sample index using:

```text
delay = total_distance / speed_of_sound

sample_index = delay × sampling_frequency
```

The sample index is then used to access the corresponding receive-channel data.

---

## 3. Phase Compensation

The delayed complex receive signal is phase compensated according to the ultrasonic carrier frequency.

The current implementation uses a 40-kHz carrier.

The phase is represented as a normalized value corresponding to one complete turn:

```text
0 → 0°
0.25 → 90°
0.50 → 180°
0.75 → 270°
1.00 → 360°
```

Sine and cosine values are generated using a lookup-table-based approximation.

Quadrant symmetry is used to reduce the required lookup-table size.

---

## 4. Beamforming

The 64 receive channels are processed in parallel.

For each image point, phase-compensated channel data are summed to generate the beamformed result.

Conceptually:

```text
Channel 1  ─┐
Channel 2  ─┤
Channel 3  ─┤
   ...       ├──► Channel Summation ───► Beamformed Pixel
Channel 64 ─┘
```

A coherence factor (CF) is also calculated from the receive-channel signals.

The coherence factor is used to reduce incoherent components and improve the quality of the beamformed image.

---

## 5. Depth-Sparse Synthetic Aperture Accumulation

A conventional synthetic aperture architecture may require storing a large partial image for every transmitter.

To reduce on-chip memory usage, this implementation uses a depth-sparse representation.

For each Tx and scan direction, only a limited number of significant depth samples are retained.

The current configuration is:

```cpp
#define K 40
```

Therefore, up to 40 partial sums are retained for each Tx and scan direction.

The partial image storage is conceptually organized as:

```text
Tx
 │
 ├── Scan Direction 1
 │      ├── Depth Candidate 1
 │      ├── Depth Candidate 2
 │      ├── ...
 │      └── Depth Candidate K
 │
 ├── Scan Direction 2
 │      └── ...
 │
 └── ...
```

The processing order is:

```cpp
for each phi
    for each theta
        for each depth
```

This allows the architecture to evaluate the depth samples associated with each scan direction before updating the sparse synthetic aperture representation.

The purpose of this approach is to reduce FPGA on-chip memory usage while retaining the dominant ultrasonic reflections.

---

## Testbench Configuration

The datasets to be tested are specified in `SA_v3_tb.cpp`.

For example:

```cpp
vector<string> datasets = {
    "0_0_100"
};
```

Additional datasets can be added as:

```cpp
vector<string> datasets = {
    "0_0_50",
    "0_0_100",
    "0_0_150",
    "30_30_50",
    "30_30_100",
    "30_30_150"
};
```

---

## File Path Configuration

The current testbench uses absolute Linux file paths.

For example:

```text
/home/hung52852/SA/testdata/<dataset_name>/
```

Before running the project on another computer, these paths should be modified to match the local directory structure.

A relative-path implementation is recommended for better portability.

For example:

```text
../TEST_DATA/<dataset_name>/
```

---

## Output Files

The testbench can generate the following files.

### Beamformed Outputs

```text
beamformed_output_tx01.txt
beamformed_output_tx02.txt
beamformed_output_tx03.txt
beamformed_output_tx04.txt
beamformed_output_tx05.txt
```

### Individual Synthetic Aperture Images

```text
SA_image_tx01.txt
SA_image_tx02.txt
SA_image_tx03.txt
SA_image_tx04.txt
SA_image_tx05.txt
```

### Summed Synthetic Aperture Image

```text
SA_image_sum.txt
```

### Top-K Results

```text
SA_topk_tx01.txt
SA_topk_tx02.txt
SA_topk_tx03.txt
SA_topk_tx04.txt
SA_topk_tx05.txt
```

Additional debugging files may also be generated during software simulation.

---

## Development Environment

The design is intended for implementation and verification using:

- AMD/Xilinx Vitis HLS
- C/C++
- Xilinx arbitrary-precision data types
- FPGA-based hardware acceleration

Main HLS libraries used include:

```cpp
#include <ap_fixed.h>
#include <ap_int.h>
#include <ap_axi_sdata.h>
#include <hls_stream.h>
#include <hls_math.h>
```

---

## Notes

This repository contains research code developed as part of a master's thesis project.

The implementation is primarily intended for:

- FPGA architecture evaluation
- Vitis HLS simulation
- Hardware synthesis
- Synthetic aperture ultrasound imaging experiments

The code is research-oriented and is not intended to be a general-purpose software beamforming library.
