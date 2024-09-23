-*- coding: utf-8 -*-
import iFeatureOmegaCLI
set = 'dairyc'
prot = iFeatureOmegaCLI.iProtein('GISAID_dairyc_prot.FASTA')
prot.display_feature_types()
prot.import_parameters('protein_params.json')
prot.get_descriptor('DPC type 1')
prot.to_csv('prot/GISAID_dairyc_prot_2mer.csv', index=True, header=True)
prot.get_descriptor('PAAC')
prot.to_csv('prot/GISAID_dairyc_prot_pseaac.csv', index=True, header=True)
prot.get_descriptor('CTriad')
prot.to_csv('prot/GISAID_dairyc_prot_ctriad.csv', index=True, header=True)
prot.get_descriptor('CTDC')
prot.to_csv('prot/GISAID_dairyc_prot_ctdc.csv', index=True, header=True)
prot.get_descriptor('CTDT')
prot.to_csv('prot/GISAID_dairyc_prot_ctdt.csv', index=True, header=True)
prot.get_descriptor('CTDD')
prot.to_csv('prot/GISAID_dairyc_prot_ctdd.csv', index=True, header=True)