# Docking Assignment

For Problem Set 6, students will perform ligand docking in an important drug target. Students are encouraged to complete the [docking tutorial](01_docking_tutorial.md) prior to starting the assignment. Students are welcome to choose any PDB entry as their starting point, but several [recommendations](#recommended-complex-structures) are listed below.

## Assignment Deliverables

1. Screenshot of spheres that represent the receptor binding site (**5 pts**)
2. Screenshot of original ligand in its original pose and a docked pose using flexible docking (**5 pts**)
3. Report on the top five ligands arising from a virtual screen of a library with at least 100 compounds (**10 pts**)
   - Screenshots of their docked poses
   - List of Descriptor scores for each top ligand
   - Rationale/discussion of which ligand the student thinks is the best starting point for drug development

## Recommended complex structures

- [7KK3](https://www.rcsb.org/structure/7KK3) - PARP1 in complex with talazoparib
- [6UOC](https://www.rcsb.org/structure/6uoc) - HDAC6 in complex with givinostat
- [8BXH](https://www.rcsb.org/structure/8bxh) - JAK2 in complex with momelotinib
- [4FS3](https://www.rcsb.org/structure/4FS3) - *S. aureus* FabI in complex with Debio-1452
- [2W9G](https://www.rcsb.org/structure/2W9G) - *S. aureus* DHFR in complex with trimethoprim
- [5YVE](https://www.rcsb.org/structure/5YVE) - P2X3 in complex with gefapixant
- [4GV1](https://www.rcsb.org/structure/4gv1) - AKT1 in complex with capivasertib
- [3S91](https://www.rcsb.org/structure/3S91) - Bromodomain in complex with JQ1
- [1NDV](https://www.rcsb.org/structure/1ndv) - Adenosine deaminase
- [1Y0X](https://www.rcsb.org/structure/1Y0X) - Thyroxine-Thyroid Hormone Receptor
- [4S0V](https://www.rcsb.org/structure/4S0V) - Orexin receptor in complex with suvorexant
- [3WZE](https://www.rcsb.org/structure/3WZE) - VEGFR in complex with sorafenib
- [2P16](https://www.rcsb.org/structure/2P16) - Factor Xa in complex with apixaban
- [6ME2](https://www.rcsb.org/structure/6ME2) - Melatonin receptor MT1 in complex with ramelteon
- [4ZUD](https://www.rcsb.org/structure/4ZUD) - Angeotensin receptor in complex with olmesartan
- [1X70](https://www.rcsb.org/structure/1X70) - Dipeptidyl peptidase IV in complex with sitagliptin
- [1HW9](https://www.rcsb.org/structure/1HW9) - HMG-CoA reductase in complex with simvastatin
- [2ZD1](https://www.rcsb.org/structure/2ZD1) - HIV reverse transcriptase in complex with rilpivirine
- [1S19](https://www.rcsb.org/structure/1S19) - Vitamin D3 receptor in complex with calcipotriol
- [1EFY](https://www.rcsb.org/structure/1EFY) - PARP in complex with benzimidazole
- [2GQG](https://www.rcsb.org/structure/2GQG) - ABL in complex with dasatinib
- [4F4P](https://www.rcsb.org/structure/4F4P) - SYK kinase
- [6UZW](https://www.rcsb.org/structure/6UZW) - NMDA glutamate receptor in complex with UBP791

## Compound Libraries

Several compound libraries have been prepared for docking and can be found in `/class/bsdrown/data/libraries`.

- `drugbank` - structures from [drugbank.ca](https://www.drugbank.ca/releases/latest#structures)
- `fragments` - collection of fragments from various vendors (e.g. Asinex, Life Chemicals, ChemDiv, Pharmablock)
- `chembridge_dvs_cl` - diversity library from the ChemBridge CORE library [ChemBridge DIVERSet-CL](https://chembridge.com/wp-content/uploads/2022/08/ChemBridge-DIVERSet-Libraries.pdf)
- `MCE_Kinase` - kinase inhibitors sold by [MedChemExpress](https://www.medchemexpress.com/)
- `LC_Kinase_Focused` - diverse set of compounds with reported kinase inhibitory activity from [Life Chemicals](https://lifechemicals.com/screening-libraries/targeted-and-focused-screening-libraries/kinase-general-libraries/kinase-focused-library)
- `LOPAC_1280` - LOPAC library sold by [SigmaAldrich](https://www.sigmaaldrich.com/US/en/product/sigma/lo1280)
