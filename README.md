# OIK-PAS-V1: Bare-Metal FPGA Spatial Voxel Screening & Allosteric TYK2 JH2 Discovery

Open-source pipeline, bare-metal Xilinx UltraScale+ PCIe BAR0 telemetry, crystallographic validation controls, and audited dual-target consensus docking artifacts contrasting **TYK2 JH2 (`PDB: 6NZP`)** against off-target **JAK2 JH1 (`PDB: 4IVA`)**.

- **Live Deployment & AWS Private Offers:** [https://oikonomia-platform.com](https://oikonomia-platform.com)
- **Preprint PDF:** [`Oikonomia_TYK2_JH2_Preprint.pdf`](./Oikonomia_TYK2_JH2_Preprint.pdf)
- **AWS EC2 F2 Golden AMI (`us-east-1`):** `ami-0a2c8c098324a5565` (`f2.6xlarge`, AFI `agfi-0b55f06c15a385187`)

---

## Hardware & Validation Summary (`agfi-0b55f06c15a385187`)
- **Silicon Substrate:** Xilinx UltraScale+ (`0000:34:00.0`) on AWS EC2 `f2.6xlarge`
- **Power & State Footprint:** `6 W` steady-state `Vccint` core power | `0.25 KB` constant `O(1)` on-chip state
- **BAR0 MMIO Latency:** `87.03 µs` P50 / `103.00 µs` P99 native C stream (`385 × 4,096-byte` `16×16×16` 3D voxel tensors)
- **Output Determinism:** `100%` bit-exact register readback (`Determinism=True`)
- **Structural Alignment & Controls:** `1.887 Å` binding-site Cα superposition RMSD (`4IVA` onto `6NZP`); `0.867` ROC-AUC separating known JH2 binders from clinical decoys

---

## Quickstart (AWS EC2 `f2.6xlarge` in `us-east-1`)

```bash
# 1. Verify live Xilinx UltraScale+ AFI & 6 W Vccint power rail
sudo fpga-describe-local-image -S 0 -R -M

# 2. Clone repository and compile userspace PCIe BAR0 streamer
git clone [https://github.com/joeiko3424-svg/oikonomia-tyk2-jh2-fpga.git](https://github.com/joeiko3424-svg/oikonomia-tyk2-jh2-fpga.git)
cd oikonomia-tyk2-jh2-fpga
gcc -O3 -march=native -shared -fPIC oikonomia_voxel_dma.c -o liboikonomia_voxel_dma.so

# 3. Run pinned dual-target screen & hardware telemetry verification
taskset -c 2 python3 run_tyk2_jh2_biorxiv.py
