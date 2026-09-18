# Thesis and Journal Paper

This repository includes two research documents related to the development of a real-time 3-D air-coupled ultrasonic imaging system based on synthetic aperture beamforming and FPGA acceleration.

The two documents share the same core research topic, but serve different purposes:

- The **master's thesis** provides the complete research background, system design, implementation details, experiments, and performance analysis.
- The **journal paper** presents a more compact and focused version of the work, emphasizing the proposed fully on-chip architecture and depth-sparse synthetic aperture reconstruction.

---

## Master's Thesis

### Depth-Sparse 3-D Synthetic Aperture Ultrasonic Imaging with a Fully On-Chip Beamforming Architecture

**Author:** Yu Xiang Hung  
**Advisor:** Dr. Bo Cheng Lai  
**Institute:** Institute of Electronics, National Yang Ming Chiao Tung University  
**Degree:** Master of Engineering in Computer Architecture  
**Date:** April 2026

The master's thesis presents the complete development of a 3-D synthetic aperture air-coupled ultrasonic imaging system implemented on an AMD/Xilinx ZCU102 FPGA platform.

The work targets real-time environment sensing applications such as autonomous vehicles, drones, smart factories, and smart warehouses. To improve the angular resolution of conventional ultrasonic sensing systems, multiple transmit events are coherently combined using synthetic aperture imaging.

The proposed hardware architecture consists of three major components:

1. **Hybrid Delay Calculator**
   - Exact transmit-distance calculation
   - Linear approximation for receive delay
   - Lookup-table-based implementation for reducing hardware complexity

2. **Streaming DAS/CF Beamformer**
   - Parallel processing across 64 receive channels
   - Quadrature demodulation
   - Delay-and-sum beamforming
   - Coherence Factor weighting
   - Fully pipelined streaming architecture

3. **Depth-Sparse Synthetic Aperture Accumulation**
   - Exploits sparsity along the depth direction
   - Retains only significant voxel responses
   - Reduces intermediate image storage
   - Enables synthetic aperture accumulation to remain fully on chip

The thesis also provides detailed discussions on system architecture, channel-level parallelism, sparse memory management, synthetic aperture scheduling, hardware implementation, Field II simulation, imaging performance, FPGA resource utilization, and throughput analysis.

The implemented system uses:

- 5 transmit events
- 64 receive channels
- 69 × 69 angular scanlines
- 284 depth samples
- ±60° field of view
- 0.3–2.0 m sensing range
- 40 kHz air-coupled ultrasound

The FPGA implementation achieves real-time 3-D ultrasonic imaging at approximately **68.56 frames/s**.

This thesis serves as the full technical documentation of the project and contains more detailed derivations, architecture descriptions, implementation methodology, and evaluation results than the journal paper.

---

## Journal Paper

### Depth-Sparse 3-D Synthetic Aperture Ultrasonic Imaging with a Fully On-Chip Beamforming Architecture

**Authors:** Yu-Xiang Hung, Po-Jui Su, Geng-Shi Jeng, and Bo-Cheng Lai  
**Target Journal:** IEEE Sensors Journal

The journal paper presents a condensed version of the system with a stronger focus on the architectural contributions required to implement 3-D synthetic aperture ultrasonic imaging entirely on FPGA on-chip memory.

The main challenge addressed by this work is the large computational and memory cost of 3-D synthetic aperture reconstruction.

Conventional synthetic aperture imaging requires storing and accumulating dense partial 3-D images from multiple transmit events. This becomes impractical for a fully on-chip FPGA implementation because of the limited amount of Block RAM available.

To solve this problem, the paper proposes a **depth-sparse synthetic aperture accumulation architecture**.

Instead of storing the complete depth profile of each scanline, only the most significant voxel responses are retained using a Top-K strategy.

For the implemented configuration:

```text
Depth samples per scanline: 284
Retained samples:           K = 30
Retained percentage:        10.56%
