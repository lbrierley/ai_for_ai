# -*- coding: utf-8 -*-

import os
import glob
import iFeatureOmegaCLI

# folders
fasta_dir = r"/hpscol02/tenant1/zoonosis-risk-ai/zoonosis-risk-ai-modelling/Protein_features_RS_2026_04_16/fastas"
out_dir = r"/hpscol02/tenant1/zoonosis-risk-ai/zoonosis-risk-ai-modelling/Protein_features_RS_2026_04_16/Results"
params_file = r"/hpscol02/tenant1/zoonosis-risk-ai/zoonosis-risk-ai-modelling/Protein_features_RS_2026_04_16/protein_params.json"

os.makedirs(out_dir, exist_ok=True)

# all fasta files in the folder
fasta_files = sorted(glob.glob(os.path.join(fasta_dir, "*.fasta")))

# descriptors to run: (descriptor name in iFeatureOmega, output suffix)
descriptors = [
    ("DPC type 1", "_2mer.csv"),
    ("PAAC", "_pseaac.csv"),
    ("CTriad", "_ctriad.csv"),
    ("CTDC", "_ctdc.csv"),
    ("CTDT", "_ctdt.csv"),
    ("CTDD", "_ctdd.csv"),
]

for fasta_file in fasta_files:
    base_name = os.path.splitext(os.path.basename(fasta_file))[0]

    print(f"Processing {base_name}...")

    prot = iFeatureOmegaCLI.iProtein(fasta_file)
    prot.import_parameters(params_file)

    for desc_name, suffix in descriptors:
        prot.get_descriptor(desc_name)
        out_file = os.path.join(out_dir, base_name + suffix)
        prot.to_csv(out_file, index=True, header=True)

    print(f"Done: {base_name}")

print("All files processed.")
