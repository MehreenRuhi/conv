# conv
Conv.py will find convergent sites between the specified species (or groups of species) given a genome sequence alignement containing all species. 
Example data are included in the "testdata" directory and can be run as follows:

for i in `ls testdata`; do python conv.py testdata/$i 'micMur1_Lepilemur_Propithecus' 'monDom5_macEug1_ponAbe2' ; done

In the above example micMur1 and monDon5 are the species for which convergent sites will be found. Lepilemur and Propithecus are the outgroup species for micMur1, and macEug1 and ponAbe2 are the outgroup species for monDom5. The names in this command must match the names in the genome sequence alignment files. Convergent sites for groups of species can be found as follows:

for i in `ls testdata`; do python conv.py testdata/$i 'micMur1_Lepilemur_Propithecus' '(colobine_rheMac2_papHam1)panTro2_gorGor1' ; done

Here the convergent sites between micMur1 and the "ancestor" of colobine, rheMac2, and papHam1 will be found. The species specified by panTro2 and gorGor1 are the outgroup for the group colobine_rheMac2_papHam. 

Two files will be created by this program, one outputting all possible convergence events between the specified groups (beginning with "whole"). The other (beginning with "conv_") will be created only if there is a convergence event and will specify the number of missing sites, the number of possible convergent sites, and the converged sites (specified by the animo acid output being given in the sequence). 
