#!/usr/bin/env python
# encoding: utf-8

# This script summarizes 95% HPD and means from a collection of trees from TreeAnnotator
# for certain nodes.
# Author: Carl H. Oliveros

import dendropy

# name of output file
output_fname = 'age.data.txt'

# List here filenames of trees annotated using TreeAnnotator 
tree_fname_list = ['rs-0-root-prior.ann.tre',
	'rs-0-root.resamp10k.ann.tre',
	'rs-0-part-prior.ann.tre',
	'rs-0-part.resamp10k.ann.tre',
	'rs-0.prior.ann.tre',
	'rs-0.resamp20k.ann.tre',
	'rs-1.resamp20k.ann.tre',
	'rs-2.resamp20k.ann.tre',
	'rs-3.resamp20k.ann.tre',
	'rs-5.resamp20k.ann.tre',
	'rs-7.resamp20k.ann.tre',
	'rs-8.resamp20k.ann.tre',
	'rs-9.resamp20k.ann.tre']

# Python dictionary of nodes of interest, defined by any two tips in the tree whose MRCA
# is the node of interest. Note that dendropy by default converts underscores to spaces
# so underscores in tip names should be represented here with spaces.
node_list = {'Acanthisittidae': ['acanthisitta chloris','xenicus gilviventris or0299521c'],
	'Creepers-sister': ['regulus regulus 6778','sitta europea 30442'],
	'Daph-Mohoua': ['daphoenositta chrysoptera 23086','mohoua albicilla 96717'],
	'Falcons': ['microhierax erythrogenys 30969','falco peregrinus'],
	'Meliph-Parda': ['meliphaga montana 12276','pardalotus striatus 8886'],
	'Menura-Atrichornis': ['menura novaehollandiae 76638','atrichornis rufescens 554508'],
	'Orthonyx-Pomatostomus': ['orthonyx temminckii 76694','pomatostomus superciliosus 8792'],																																					
	'Passerines-Parrots': ['passer domesticus 28860','psittacus erithacus'],	
	'Peltops': ['peltops blainvillii 809112','strepera graculina 9660'],
	'Strigops-Nestor': ['strigops habroptila 23361','nestor notabilis'],
	'Osc-Suboscines': ['passer domesticus 28860','tyrannus albogularis b7247'],
	'Root': ['cariama cristata','passer domesticus 28860'],
	'OW-NW Suboscines': ['eurylaimus ochromalus b50329','tyrannus albogularis b7247'],
	'Crown Oscines':['menura novaehollandiae 76638','passer domesticus 28860'],
	'Corv-Pass-split': ['corvus corax 30042','passer domesticus 28860'],
	'Bombycilloidea': ['bombycilla garrulus 21744','dulus dominicus 6414'],
	'Tyrranida': ['tyrannus albogularis b7247','pipra filicauda b4330'],
	'Corvoidea': ['corvus corax 30042','rhipidura javanica 17717'],
	'Sylviida': ['curruca nana 28718','hyliota flavigaster 439359'],
	'Core-Passeridans': ['peucedramus taeniatus 9396','passer domesticus 28860']}
	
outfile = open(output_fname, 'w')
outfile.write('node\tfilename\t95high\t95low\tmean\n')

for tree_fname in tree_fname_list:
	tree = dendropy.Tree.get_from_path(tree_fname,'nexus')
	for nd in node_list.keys():
		node = tree.mrca(taxon_labels = node_list[nd])
		height = ''
		low = ''
		high = ''
		for item in node._annotations:
			if item.name == 'height':
				height = item.value
			elif item.name == 'height_95%_HPD':
				low = item.value[0]
				high = item.value[1]
		outfile.write('{0}\t{1}\t{2}\t{3}\t{4}\n'.format(nd, tree_fname, high, low, height))
		
outfile.close()


