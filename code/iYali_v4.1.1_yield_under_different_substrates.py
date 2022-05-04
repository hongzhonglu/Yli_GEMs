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
cobra_config.solver = "glpk_exact"
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
result_df["yield(gDW/mmol C)"] = result_df["growth_rate"]/18
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

# only acetate
biomass_acetate = []
for x, y in zip(for_flux, ac_flux):
    model3.reactions.get_by_id("y001714").bounds = (0, 0)  # D-glucose exchange
    model3.reactions.get_by_id("y001793").bounds = (0, 0)  # formate exchange
    model3.reactions.get_by_id("y001634").bounds = (y, 0)  # acetate exchange
    solution10 = model3.optimize()
    biomass_acetate_v = solution10.objective_value
    biomass_acetate.append(biomass_acetate_v)

result_df2 = pd.DataFrame({"formate_uptake":for_flux,"acetate_uptake":ac_flux,"growth_rate":biomass, "growth_rate_acetate":biomass_acetate})
result_df2["yield(gDW/mmol C)"] = result_df2["growth_rate"]/18
result_df2.to_excel("result/formate_acetate_for_yli.xlsx")


# add
model1 = model.copy()

# L-tyrosine[c] => ammonium[c] + trans-4-coumarate[c]
# trans-4-coumarate[c] + ATP[c] + coenzyme A[c] => AMP[c] + diphosphate[c] + 4-coumaroyl-CoA[c]
# 3 H+[c] + 4-coumaroyl-CoA[c] + 3 malonyl-CoA[c] => 4 coenzyme A[c] + 4 carbon dioxide[c] + resveratrol[c]
# resveratrol[c] => resveratrol[e]
# resveratrol[e] =>


from cobra import Reaction, Metabolite


# firstly add the new metabolite
model1.add_metabolites(Metabolite('n_0001', compartment='c', formula='C9H7O3', name='trans-4-coumarate'))
model1.add_metabolites(Metabolite('n_0002', compartment='c', formula='C30H38N7O18P3S', name='4-coumaroyl-CoA'))
model1.add_metabolites(Metabolite('n_0003', compartment='c', formula='C14H12O3', name='resveratrol'))
model1.add_metabolites(Metabolite('n_0004', compartment='e', formula='C14H12O3', name='resveratrol'))



# reaction 1
dict1 = {model1.metabolites.get_by_id('s_1051'): -1,
         model1.metabolites.get_by_id('s_0419'):1,
         model1.metabolites.get_by_id('n_0001'):1
        }

dict2 = {model1.metabolites.get_by_id('n_0001'): -1,
         model1.metabolites.get_by_id('s_0434'):-1,
         model1.metabolites.get_by_id('s_0529'):-1,
         model1.metabolites.get_by_id('s_0423'):1,
         model1.metabolites.get_by_id('s_0633'):1,
         model1.metabolites.get_by_id('n_0002'):1
        }

dict3 = {model1.metabolites.get_by_id('s_0794'): -3,
         model1.metabolites.get_by_id('n_0002'):-1,
         model1.metabolites.get_by_id('s_1101'):-3,
         model1.metabolites.get_by_id('s_0529'):4,
         model1.metabolites.get_by_id('s_0456'):4,
         model1.metabolites.get_by_id('n_0003'):1
        }

dict4 = {model1.metabolites.get_by_id('n_0003'): -1,
         model1.metabolites.get_by_id('n_0004'):1
        }


# exchange
dict5 = {model1.metabolites.get_by_id('n_0004'): -1
        }



# using the general procedures to add the new function
reaction = Reaction('new1')
reaction.name = 'new1'
reaction.subsystem = 'new_added'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.EC = ''
reaction.add_metabolites(dict1)
reaction.gene_reaction_rule = ''
model1.add_reactions([reaction])

reaction = Reaction('new2')
reaction.name = 'new2'
reaction.subsystem = 'new_added'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.EC = ''
reaction.add_metabolites(dict2)
reaction.gene_reaction_rule = ''
model1.add_reactions([reaction])


reaction = Reaction('new3')
reaction.name = 'new3'
reaction.subsystem = 'new_added'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.EC = ''
reaction.add_metabolites(dict3)
reaction.gene_reaction_rule = ''
model1.add_reactions([reaction])

reaction = Reaction('new4')
reaction.name = 'new4'
reaction.subsystem = 'new_added'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.EC = ''
reaction.add_metabolites(dict4)
reaction.gene_reaction_rule = ''
model1.add_reactions([reaction])

reaction = Reaction('new5')
reaction.name = 'new5'
reaction.subsystem = 'new_added'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.EC = ''
reaction.add_metabolites(dict5)
reaction.gene_reaction_rule = ''
model1.add_reactions([reaction])

for rxn in model1.reactions:
    print(rxn.id)


# simulation
# first fix the growth rate
model1.objective = 'new5'
solution = model1.optimize()
product_v = solution.objective_value
model1.reactions.get_by_id("y002111").bounds = (0.1, 0.1)
solution = model1.optimize()
product_v = solution.objective_value






# how to combine
for_flux = np.linspace(-18, 0, 30)
ac_flux = [(-18-x)/2 for x in for_flux]
product = []
for x, y in zip(for_flux, ac_flux):
    model1.reactions.get_by_id("y001714").bounds = (0, 0)  # D-glucose exchange
    model1.reactions.get_by_id("y001793").bounds = (x, 0)  # formate exchange
    model1.reactions.get_by_id("y001634").bounds = (y, 0)  # acetate exchange
    solution10 = model1.optimize()
    product_v = solution10.objective_value
    product.append(product_v)

# only acetate
product_acetate = []
for x, y in zip(for_flux, ac_flux):
    model1.reactions.get_by_id("y001714").bounds = (0, 0)  # D-glucose exchange
    model1.reactions.get_by_id("y001793").bounds = (0, 0)  # formate exchange
    model1.reactions.get_by_id("y001634").bounds = (y, 0)  # acetate exchange
    solution10 = model1.optimize()
    product_acetate_v = solution10.objective_value
    product_acetate.append(product_acetate_v)

result_df2 = pd.DataFrame({"formate_uptake":for_flux,"acetate_uptake":ac_flux,"product_rate":product, "product_rate_acetate":product_acetate})
result_df2["product_yield(gDW/mmol C)"] = result_df2["product_rate"]/18
result_df2["product_yield_acetate(gDW/mmol C)"] = result_df2["product_rate_acetate"]/(-2*result_df2["acetate_uptake"])

result_df2.to_excel("result/product_formate_acetate_for_yli.xlsx")










# glucose exchange x acetate exchange
model1.add_boundary(model.metabolites.get_by_id("m1641"), type="sink")
model1.objective = 'SK_m1641'
model1.objective = 'new5'
model1.reactions.get_by_id("y002111").bounds = (0, 0)
model1.reactions.get_by_id("y001714").bounds = (-0.5, 0) # D-glucose exchange
model1.reactions.get_by_id("y001793").bounds = (0, 0) # formate exchange
model1.reactions.get_by_id("y001634").bounds = (-1.5, 0) # acetate exchange
solution10 = model1.optimize()
# how to combine
gluc_flux = np.linspace(-3, 0, 30)
ac_flux = [(-3-x)*3 for x in gluc_flux]
product = []
for x, y in zip(gluc_flux, ac_flux):
    model1.reactions.get_by_id("y001714").bounds = (x, 0)  # D-glucose exchange
    model1.reactions.get_by_id("y001793").bounds = (0, 0)  # formate exchange
    model1.reactions.get_by_id("y001634").bounds = (y, 0)  # acetate exchange
    solution10 = model1.optimize()
    product_v = solution10.objective_value
    product.append(product_v)

result_df = pd.DataFrame({"glucose_uptake":gluc_flux,"acetate_uptake":ac_flux,"product_rate":product})
result_df["yield(gDW/mmol C)"] = result_df["product_rate"]/18
result_df.to_excel("result/product_glucose_acetate_for_yli.xlsx")








