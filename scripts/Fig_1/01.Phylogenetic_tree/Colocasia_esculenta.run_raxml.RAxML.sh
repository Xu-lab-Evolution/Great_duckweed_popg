#! /bin/bash

raxml-ng --msa ../../../data//Fig_1/01.Phylogenetic_tree/Colocasia_esculenta.run_raxml/Colocasia_esculenta.run_raxml.fa --model TVM+G4  --all  --threads 40 --bs-metric fbp,tbe  --seed 142857 --tree pars{50},rand{50} --prefix ../../../data//Fig_1/01.Phylogenetic_tree/Colocasia_esculenta.run_raxml/Colocasia_esculenta.run_raxml.raxml
