# -*- coding: utf-8 -*-
# -*- python 3 -*-
# -*- hongzhong Lu -*-

# to do
# update biomass composition for 4 strains       reaction ID: xBIOMASS
# update GEMs based on metabolic engineering strategy
# simulate the phenotype and get flux data





from cobra.io import read_sbml_model, load_matlab_model, write_sbml_model
import cobra
cobra_config = cobra.Configuration()


model = read_sbml_model("model/iYali_v4.1.1.xml")
rxn_inf= []
# print rxn
for rxn in model.reactions:
    rxn_inf.append(rxn.name)
    if rxn.name == "bicarbonate exchange":
        print(rxn.id)



# print gene
for gene in model.genes:
    print(gene.id)


exchange_rxn = [x for x in rxn_inf if "exchange" in x]
# D-glucose exchange; bicarbonate exchange; oxygen exchange

# test using the new data from Fujing
model = read_sbml_model("model/iYali_v4.1.1.xml")

# strain 1
model.reactions.get_by_id("y001714").bounds = (-1.04, 0) # D-glucose exchange
#model.reactions.get_by_id("R_EX_o2(e)").bounds = (-2,0) # current no data
model.reactions.get_by_id("y001663").bounds = (2.06, 2.06) # bicarbonate exchange
solution1 = model.optimize()

# strain 2
model.reactions.get_by_id("y001714").bounds = (-1.34, 0)
#model.reactions.get_by_id("R_EX_o2(e)").bounds = (-2,0) # current no data
model.reactions.get_by_id("y001663").bounds = (2.29, 2.29)
solution2 = model.optimize()

# strain 3
model.reactions.get_by_id("y001714").bounds = (-1.5, 0)
#model.reactions.get_by_id("R_EX_o2(e)").bounds = (-2,0) # current no data
model.reactions.get_by_id("y001663").bounds = (3.13, 3.13)
solution3 = model.optimize()

# strain 4
model.reactions.get_by_id("y001714").bounds = (-1.237, 0)
#model.reactions.get_by_id("R_EX_o2(e)").bounds = (-2,0) # current no data
model.reactions.get_by_id("y001663").bounds = (3.00, 3.00)
solution4 = model.optimize()

# strain 5
model.reactions.get_by_id("y001714").bounds = (-1.15, 0)
#model.reactions.get_by_id("R_EX_o2(e)").bounds = (-2,0) # current no data
model.reactions.get_by_id("y001663").bounds = (2.93, 2.93)
solution5 = model.optimize()



# test using Ed data
model = read_sbml_model("model/iYali_v4.1.1.xml")
model.reactions.get_by_id("y001714").bounds = (-0.61, 0)
#model.reactions.get_by_id("R_EX_o2(e)").bounds = (-2,0) # current no data
model.reactions.get_by_id("y001663").bounds = (1.5, 1.5)
solution = model.optimize()