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
import pandas as pd
import numpy as np
model = read_sbml_model("model/iYali_v4.1.1.xml")


id_all = []
name_all = []
formula_all = []
upper_all = []
lower_all = []
for rxn in model.reactions:
    print(rxn.id)
    id_all.append(rxn.id)
    name_all.append(rxn.name)
    formula_all.append(rxn.reaction)
    upper_all.append(rxn.upper_bound)
    lower_all.append(rxn.lower_bound)

rxn_df = pd.DataFrame({"id":id_all,"name":name_all,"formula":formula_all,"lower_bound":lower_all, "upper_bound":upper_all})
rxn_df.to_excel("result/iYali_v4.1.1.xlsx")
rxn_df_exchange = rxn_df[rxn_df["name"].str.contains("exchange")]

# y001793 formate exchange
model1 = model.copy()
model1.reactions.get_by_id("y001714").bounds = (0, 0) # D-glucose exchange
model1.reactions.get_by_id("y001793").bounds = (-6, 0) # formate exchange
solution2 = model1.optimize()

# y001634 acetate exchange
model2 = model.copy()
model2.reactions.get_by_id("y001714").bounds = (0, 0) # D-glucose exchange
model2.reactions.get_by_id("y001793").bounds = (0, 0) # formate exchange
model2.reactions.get_by_id("y001634").bounds = (-3, 0) # acetate exchange
solution3 = model2.optimize()

# co-utilization
modelx = model.copy()
modelx.reactions.get_by_id("y001714").bounds = (0, 0) # D-glucose exchange
modelx.reactions.get_by_id("y001793").bounds = (-3, 0) # formate exchange
modelx.reactions.get_by_id("y001634").bounds = (-1.5, 0) # acetate exchange
solution4 = modelx.optimize()


# glucose exchange
model3 = model.copy()
model3.reactions.get_by_id("y001714").bounds = (-1, 0) # D-glucose exchange
model3.reactions.get_by_id("y001793").bounds = (0, 0) # formate exchange
model3.reactions.get_by_id("y001634").bounds = (0, 0) # acetate exchange
solution1 = model3.optimize()


# glucose exchange x acetate exchange
model3 = model.copy()
model3.reactions.get_by_id("y001714").bounds = (-0.5, 0) # D-glucose exchange
model3.reactions.get_by_id("y001793").bounds = (0, 0) # formate exchange
model3.reactions.get_by_id("y001634").bounds = (-1.5, 0) # acetate exchange
solution10 = model3.optimize()
# how to combine
gluc_flux = np.linspace(-3, 0, 30)
ac_flux = [(-3-x)*3 for x in gluc_flux]
biomass = []
for x, y in zip(gluc_flux, ac_flux):
    model3.reactions.get_by_id("y001714").bounds = (x, 0)  # D-glucose exchange
    model3.reactions.get_by_id("y001793").bounds = (0, 0)  # formate exchange
    model3.reactions.get_by_id("y001634").bounds = (y, 0)  # acetate exchange
    solution10 = model3.optimize()
    biomass_v = solution10.objective_value
    biomass.append(biomass_v)

result_df = pd.DataFrame({"glucose_uptake":gluc_flux,"acetate_uptake":ac_flux,"growth_rate":biomass})
result_df.to_excel("result/glucose_acetate_for_yli.xlsx")



# how to combine
for_flux = np.linspace(-18, 0, 30)
ac_flux = [(-18-x)/2 for x in for_flux]
biomass = []
for x, y in zip(for_flux, ac_flux):
    model3.reactions.get_by_id("y001714").bounds = (0, 0)  # D-glucose exchange
    model3.reactions.get_by_id("y001793").bounds = (x, 0)  # formate exchange
    model3.reactions.get_by_id("y001634").bounds = (y, 0)  # acetate exchange
    solution10 = model3.optimize()
    biomass_v = solution10.objective_value
    biomass.append(biomass_v)

result_df2 = pd.DataFrame({"formate_uptake":for_flux,"acetate_uptake":ac_flux,"growth_rate":biomass})
result_df2.to_excel("result/formate_acetate_for_yli.xlsx")












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