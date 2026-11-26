# OIK-PAS-V1: Direct Userspace PCIe BAR0 Spatial Voxel Screening & Allosteric TYK2 JH2 Discovery
Open-source hardware-gated deterministic pipeline, direct userspace Xilinx UltraScale+ PCIe BAR0 MMIO telemetry (kernel-bypass on AWS EC2 f2.6xlarge), crystallographic validation controls, and audited dual-target consensus docking artifacts contrasting allosteric **TYK2 JH2 (PDB: 6NZP)** against catalytic **JAK2 JH1 (PDB: 4IVA)**.

* **Live Deployment & Qualification Sprints:** [https://oikonomia-platform.com](https://oikonomia-platform.com)
* **Preprint PDF:** `Oikonomia_TYK2_JH2_ChemRxiv.pdf`

## Hardware & Validation Summary
* **Silicon Substrate:** Xilinx UltraScale+ (`0000:34:00.0`) on single-tenant AWS EC2 `f2.6xlarge` (accessed via direct userspace PCIe BAR0 MMIO)
* **Power & State Footprint:** 6 W steady-state $V_{CCINT}$ core power | 0.25 KB constant $O(1)$ on-chip state
* **BAR0 MMIO Latency:** 87.03 µs P50 / 103.00 µs P99 native C stream (385 × 4,096-byte `16×16×16` 3D voxel tensors)
* **Output Determinism:** 100% bit-exact register readback (`Determinism=True`)
* **Structural Alignment & Controls:** 1.887 Å binding-site Cα superposition RMSD (`4IVA` onto `6NZP`); 0.867 ROC-AUC separating known JH2 binders from clinical decoys

> **Note on Public Benchmark Scope:** This repository characterizes direct userspace PCIe BAR0 MMIO latency (~87.03 µs / 4 KB at 6 W core power, 0.25 KB state) and deterministic state verification alongside an open-source CPU dual-target Vina docking workflow yielding the selective allosteric lead candidate **OIK-TYK2-003**. Proprietary RTL scoring pipelines remain encapsulated inside Oikonomia's encrypted hardware enclave.

## Key Artifacts
* `run_tyk2_jh2_biorxiv.py` & `oikonomia_voxel_dma.c` — Sequence-aligned receptor superposition, RDKit/Meeko preparation, bare-metal C BAR0 voxel streamer, and 3-seed AutoDock Vina consensus runner.
* `top_tyk2_jh2_hits_audited.csv` — Full physicochemical, isoform selectivity ($\Delta\Delta G_{\text{TYK2}-\text{JAK2}}$), and hinge/pocket contact table (`Val690`, `Glu688`, `Ala671`, `Lys642`).
* `complex_TYK2_rank1_OIK-TYK2-003.pdb` & `rank2–rank6.pdb` — PyMOL-ready 3D complexes.
* `hardware_telemetry.log` & `biorxiv_tyk2_jh2_summary.json` — Raw silicon latency logs and audited summary JSON.

## Zero-Integration 72-Hour Dual-Target Selectivity Sprint
Computational chemistry and structural biology teams do not need to provision AWS F2 instances or configure PCIe BAR0 drivers. Oikonomia executes custom isoform-selectivity screens inside our isolated, single-tenant hardware vault:

1. **Submit Target Inputs:** Email two public target PDB IDs (on-target vs. off-target) and your `SMILES` library to `contact@oikonomia-platform.com`.
2. **Private Vault Execution:** Automated PAINS/Brenk triage, `87.03 µs` BAR0 spatial voxel screening (`6 W` steady-state), and multi-seed consensus docking execute inside Oikonomia's air-gapped AWS F2 enclave.
3. **Cryptographically Signed Dossier:** Receive ranked selectivity candidates ($\Delta\Delta G$, polar contact geometry, PyMOL complexes) and `SHA-256` hardware verification receipts within 72 hours (`$7,500–$12,000` flat sprint via AWS Marketplace Private Offer or invoice).

## Disclaimer
All results represent validated userspace BAR0 hardware characterization and in silico computational docking candidates. No wet-lab synthesis or biochemical assays were performed.
