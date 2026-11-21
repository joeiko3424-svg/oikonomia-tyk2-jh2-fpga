with open("/home/ubuntu/tyk2_jh2_biorxiv/run_tyk2_jh2_biorxiv.py", "r") as f:
    src = f.read()

old_align = """    fixed_ca, moving_ca = [], []
    for (s1, e1), (s2, e2) in zip(*aln.aligned):
        for i1, i2 in zip(range(s1, e1), range(s2, e2)):
            if seq_tyk2[i1] == seq_jak2[i2]:
                fixed_ca.append(tyk2_res[i1]["CA"])
                moving_ca.append(jak2_res[i2]["CA"])"""

new_align = """    fixed_ca, moving_ca = [], []
    for (s1, e1), (s2, e2) in zip(*aln.aligned):
        for i1, i2 in zip(range(s1, e1), range(s2, e2)):
            if seq_tyk2[i1] == seq_jak2[i2]:
                if np.linalg.norm(tyk2_res[i1]["CA"].get_coord() - pocket_center) <= 12.5:
                    fixed_ca.append(tyk2_res[i1]["CA"])
                    moving_ca.append(jak2_res[i2]["CA"])"""

src = src.replace('--size_x", "22", "--size_y", "22", "--size_z", "22"',
                  '--size_x", "16", "--size_y", "16", "--size_z", "16"')

old_prot_a = 'prot_a_aff, out_a = run_vina_docking(tyk2_pdbqt, pdbqt_a, "CTRL_DEUCRA_ProtA", center, exhaustiveness=16, seed=42)'
new_prot_a = """lig_in_a = os.path.join(WORK_DIR, "CTRL_DEUCRA_ProtA_s42_in.pdbqt")
            out_a = os.path.join(WORK_DIR, "CTRL_DEUCRA_ProtA_s42_out.pdbqt")
            with open(lig_in_a, "w") as fa:
                fa.write(pdbqt_a)
            res_a = subprocess.run(
                ["vina", "--receptor", tyk2_pdbqt, "--ligand", lig_in_a, "--local_only", "--out", out_a],
                capture_output=True, text=True, check=True
            )
            for ln in res_a.stdout.splitlines():
                if any(k in ln for k in ("Estimated Free Energy of Binding", "Affinity:")) or re.match(r"^\\s*1\\s+", ln):
                    nums = re.findall(r"-?\\d+\\.\\d+", ln)
                    if nums:
                        prot_a_aff = float(nums[0])
                        break"""

assert old_align in src and old_prot_a in src, "Patch anchor mismatch"
src = src.replace(old_align, new_align, 1).replace(old_prot_a, new_prot_a, 1)

with open("/home/ubuntu/tyk2_jh2_biorxiv/run_tyk2_jh2_biorxiv.py", "w") as f:
    f.write(src)
print("[+] Hardened TYK2 JH2 pocket alignment and Protocol A crystallographic control")
