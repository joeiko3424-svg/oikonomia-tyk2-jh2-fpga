# OIK-PAS-V1: Bare-Metal FPGA Spatial Voxel Screening & Allosteric TYK2 JH2 Discovery

Open-source pipeline, bare-metal Xilinx UltraScale+ PCIe BAR0 telemetry, crystallographic validation controls, and audited dual-target consensus docking artifacts (`TYK2 JH2` PDB `6NZP` vs. off-target `JAK2` PDB `4IVA`).

## Hardware & Validation Summary (`agfi-0b55f06c15a385187`)
- **Silicon Substrate:** Xilinx UltraScale+ (`0000:34:00.0`) on AWS EC2 `f2.6xlarge`
- **On-Chip State Footprint:** 0.25 KB ($O(1)$ constant) | **Output Determinism:** 100% Bit-Exact (`Determinism=True`)
- **Latency:** 87.03 µs P50 / 103.00 µs P99 native C BAR0 stream (`385 × 4 KB` voxel tensors)
- **Structural Alignment & Controls:** 1.887 Å binding-site Cα superposition RMSD (`4IVA` onto `6NZP`); **0.867 ROC-AUC** active/decoy separation

## Key Artifacts
- `run_tyk2_jh2_biorxiv.py` & `oikonomia_voxel_dma.c`: Sequence-aligned receptor superposition, RDKit/Meeko preparation, bare-metal C BAR0 voxel streamer, and 3-seed AutoDock Vina consensus runner.
- `top_tyk2_jh2_hits_audited.csv`: Full physicochemical, isoform selectivity ($\Delta\Delta G_{\text{TYK2}-\text{JAK2}}$), and hinge/pocket contact table (`Val690`, `Glu688`, `Ala671`, `Lys642`).
- `complex_TYK2_rank1_OIK-TYK2-003.pdb` & `rank2`–`rank6.pdb`: PyMOL-ready 3D complexes.

## Disclaimer
All results are *in silico* computational and hardware benchmarks. No wet-lab synthesis or biochemical assays were performed.
