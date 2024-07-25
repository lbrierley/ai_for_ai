-*- coding: utf-8 -*-

import iFeatureOmegaCLI

prot = iFeatureOmegaCLI.iProtein('/home/lbrier/S3/data/full/GISAID_human_prot.FASTA')

prot.display_feature_types()
prot.import_parameters('/home/lbrier/S3/scripts/protein_params.json')

prot.get_descriptor("DPC type 1")
prot.to_csv("/home/lbrier/S3/data/full/GISAID_human_prot_2mer.csv", index=True, header=True)

prot.get_descriptor("TPC type 1")
prot.to_csv("/home/lbrier/S3/data/full/GISAID_human_prot_3mer.csv", index=True, header=True)

prot.get_descriptor("PAAC")
prot.to_csv("/home/lbrier/S3/data/full/GISAID_human_prot_pseaac.csv", index=True, header=True)

prot.get_descriptor("CTriad")
prot.to_csv("/home/lbrier/S3/data/full/GISAID_human_prot_ctriad.csv", index=True, header=True)

prot.get_descriptor("CTDC")
prot.to_csv("/home/lbrier/S3/data/full/GISAID_human_prot_ctdc.csv", index=True, header=True)

prot.get_descriptor("CTDT")
prot.to_csv("/home/lbrier/S3/data/full/GISAID_human_prot_ctdt.csv", index=True, header=True)

prot.get_descriptor("CTDD")
prot.to_csv("/home/lbrier/S3/data/full/GISAID_human_prot_ctdd.csv", index=True, header=True)



prot = iFeatureOmegaCLI.iProtein('/home/lbrier/S3/data/full/NCBI_human_prot.FASTA')
prot.import_parameters('/home/lbrier/S3/scripts/protein_params.json')

prot.get_descriptor("DPC type 1")
prot.to_csv("/home/lbrier/S3/data/full/NCBI_human_prot_2mer.csv", index=True, header=True)

prot.get_descriptor("PAAC")
prot.to_csv("/home/lbrier/S3/data/full/NCBI_human_prot_pseaac.csv", index=True, header=True)

prot.get_descriptor("CTriad")
prot.to_csv("/home/lbrier/S3/data/full/NCBI_human_prot_ctriad.csv", index=True, header=True)

prot.get_descriptor("CTDC")
prot.to_csv("/home/lbrier/S3/data/full/NCBI_human_prot_ctdc.csv", index=True, header=True)

prot.get_descriptor("CTDT")
prot.to_csv("/home/lbrier/S3/data/full/NCBI_human_prot_ctdt.csv", index=True, header=True)

prot.get_descriptor("CTDD")
prot.to_csv("/home/lbrier/S3/data/full/NCBI_human_prot_ctdd.csv", index=True, header=True)



prot = iFeatureOmegaCLI.iProtein('/home/lbrier/S3/data/full/GISAID_avian_prot.FASTA')
prot.import_parameters('/home/lbrier/S3/scripts/protein_params.json')

prot.get_descriptor("DPC type 1")
prot.to_csv("/home/lbrier/S3/data/full/GISAID_avian_prot_2mer.csv", index=True, header=True)

prot.get_descriptor("PAAC")
prot.to_csv("/home/lbrier/S3/data/full/GISAID_avian_prot_pseaac.csv", index=True, header=True)

prot.get_descriptor("CTriad")
prot.to_csv("/home/lbrier/S3/data/full/GISAID_avian_prot_ctriad.csv", index=True, header=True)

prot.get_descriptor("CTDC")
prot.to_csv("/home/lbrier/S3/data/full/GISAID_avian_prot_ctdc.csv", index=True, header=True)

prot.get_descriptor("CTDT")
prot.to_csv("/home/lbrier/S3/data/full/GISAID_avian_prot_ctdt.csv", index=True, header=True)

prot.get_descriptor("CTDD")
prot.to_csv("/home/lbrier/S3/data/full/GISAID_avian_prot_ctdd.csv", index=True, header=True)



prot = iFeatureOmegaCLI.iProtein('/home/lbrier/S3/data/full/NCBI_avian_prot.FASTA')
prot.import_parameters('/home/lbrier/S3/scripts/protein_params.json')

prot.get_descriptor("DPC type 1")
prot.to_csv("/home/lbrier/S3/data/full/NCBI_avian_prot_2mer.csv", index=True, header=True)

prot.get_descriptor("PAAC")
prot.to_csv("/home/lbrier/S3/data/full/NCBI_avian_prot_pseaac.csv", index=True, header=True)

prot.get_descriptor("CTriad")
prot.to_csv("/home/lbrier/S3/data/full/NCBI_avian_prot_ctriad.csv", index=True, header=True)

prot.get_descriptor("CTDC")
prot.to_csv("/home/lbrier/S3/data/full/NCBI_avian_prot_ctdc.csv", index=True, header=True)

prot.get_descriptor("CTDT")
prot.to_csv("/home/lbrier/S3/data/full/NCBI_avian_prot_ctdt.csv", index=True, header=True)

prot.get_descriptor("CTDD")
prot.to_csv("/home/lbrier/S3/data/full/NCBI_avian_prot_ctdd.csv", index=True, header=True)