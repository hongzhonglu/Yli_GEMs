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
from cobra import  Reaction, Metabolite



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


# print rxn
rxn_inf= []
for rxn in model.reactions:
    rxn_inf.append(rxn.name)
    if rxn.name == "bicarbonate exchange":
        print(rxn.id)

# print gene
for gene in model.genes:
    print(gene.id)


exchange_rxn = [x for x in rxn_inf if "exchange" in x]


# First update the model
# firstly add the new metabolite
model.add_metabolites(Metabolite('s_9000', compartment='c', formula='C5H4O4', name='itaconate'))
model.add_metabolites(Metabolite('s_9001', compartment='e', formula='C5H4O4', name='itaconate'))
for met in model.metabolites:
    print(met.id)


# define the two reactions
# reaction cis-aconitate[m] => cis-aconitate[c]
dict1 = {model.metabolites.get_by_id('s_0517'): -1,
         model.metabolites.get_by_id('s_0516'): 1,
        }

# reaction cis-aconitate[c] + H[c]=>  itaconate[c] + CO2[c]
dict2 = {model.metabolites.get_by_id('s_0516'): -1,
         model.metabolites.get_by_id('s_0794'): -1,
         model.metabolites.get_by_id('s_9000'): 1,
         model.metabolites.get_by_id('s_0456'): 1,
        }


# reaction itaconate[c] => itaconate[e]
dict3 = {model.metabolites.get_by_id('s_9000'): -1,
         model.metabolites.get_by_id('s_9001'): 1,
        }

#  sink reaction
dict4 = {model.metabolites.get_by_id('s_9001'): -1
        }


# add the reactions into the model
# using the general procedures to add the new function
reaction = Reaction('new_aconitate_trasnport')
reaction.name = 'new_aconitate_trasnport'
reaction.subsystem = 'new added'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.EC = ''
reaction.add_metabolites(dict1)
reaction.gene_reaction_rule = 'AtMTT'
model.add_reactions([reaction])


# 2
# using the general procedures to add the new function
reaction = Reaction('new_cis_aconitate_decarboxylase')
reaction.name = 'new_cis_aconitate_decarboxylase'
reaction.subsystem = 'new added'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.EC = ''
reaction.add_metabolites(dict2)
reaction.gene_reaction_rule = 'AtCAD'
model.add_reactions([reaction])


# 3
reaction = Reaction('new_itaconate_trasnport')
reaction.name = 'new_itaconate_trasnport'
reaction.subsystem = 'new added'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.EC = ''
reaction.add_metabolites(dict3)
reaction.gene_reaction_rule = ''
model.add_reactions([reaction])

# 4
reaction = Reaction('exchange_itaconate')
reaction.name = 'exchange_itaconate'
reaction.subsystem = 'new added'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.EC = ''
reaction.add_metabolites(dict4)
reaction.gene_reaction_rule = ''
model.add_reactions([reaction])

# check the reaction
all_rxn = []
for rxn in model.reactions:
    print(rxn.id)
    all_rxn.append(rxn.id)

# simulation for the check
model.objective = 'exchange_itaconate'
model.reactions.get_by_id("y001714").bounds = (-1.14, 0) # D-glucose exchange
solution1 = model.optimize()





# glucose exchange
model3 = model.copy()


# remove gene - YALI0F06578g catalyze reactions including y000127 y000132
model3.reactions.get_by_id("y000127").bounds = (-1000, 1000) # reaction deletion
model3.reactions.get_by_id("y000132").bounds = (-1000, 1000) # reaction deletion


# remove gene - YALI0E16797g catalyze reactions including y102884 y102948
model3.reactions.get_by_id("y102884").bounds = (0, 0) # reaction deletion
model3.reactions.get_by_id("y102948").bounds = (0, 0) # reaction deletion

# remove gene - YALI0E32769g catalyze reactions including y000336 y000336_1
model3.reactions.get_by_id("y000336").bounds = (-1000, 1000) # reaction deletion
model3.reactions.get_by_id("y000336_1").bounds = (-1000, 1000) # reaction deletion

# remove gene - YALI0D07986g can't find it in the model


# simulation
model3.objective = 'y002111' # growth
model3.reactions.get_by_id("y001714").bounds = (-1.14, 0) # D-glucose exchange
model.reactions.get_by_id("y001663").bounds = (2.9, 2.9) # bicarbonate exchange
model.reactions.get_by_id("y001992").bounds = (-6.34, -6.34) # oxygen exchange
# by-product
model3.reactions.get_by_id("y001793").bounds = (0, 0) # formate exchange
model3.reactions.get_by_id("y001634").bounds = (0, 0) # acetate exchange
model3.reactions.get_by_id("y001687").bounds = (0.003, 0.003) # citrate exchange
model3.reactions.get_by_id("y002033").bounds = (0.0086, 0.0086) # Pyruvate exchange
model3.reactions.get_by_id("y001761").bounds = (0.0144, 0.0144) # Ethanol exchange
model3.reactions.get_by_id("y002056").bounds = (0.0115, 0.0115) # Succinate exchange
model3.reactions.get_by_id("exchange_itaconate").bounds = (0.117, 0.117) # itaconate exchange
solution1 = model3.optimize()


# import the excel file from fujing
all_data = pd.read_excel("data/fujing_dataset/Results sample_Data_20221011.xlsx")
all_data = all_data.iloc[0:15, 0:17]


all_fluxes = pd.DataFrame({"rxnID": all_rxn})

for i, x in all_data.iterrows():
    #print(i)
    #print(i, x)
    qs = x[1]
    qco2 = x[2]
    qo2 = x[3]
    qcit = x[4]
    qpyr = x[5]
    qeth = x[6]
    qsucc = x[7]
    qita = x[8]
    model3.objective = 'y002111'  # growth
    model3.reactions.get_by_id("y001714").bounds = (qs*(-1), 0)  # D-glucose exchange
    model3.reactions.get_by_id("y001663").bounds = (qco2, qco2)  # bicarbonate exchange
    model3.reactions.get_by_id("y001992").bounds = (qo2*(-1), qo2*(-1))  # oxygen exchange
    # by-product
    model3.reactions.get_by_id("y001793").bounds = (0, 0)  # formate exchange
    model3.reactions.get_by_id("y001634").bounds = (0, 0)  # acetate exchange
    model3.reactions.get_by_id("y001687").bounds = (x[4], x[4])  # citrate exchange
    model3.reactions.get_by_id("y002033").bounds = (x[5], x[5])  # Pyruvate exchange
    model3.reactions.get_by_id("y001761").bounds = (x[6], x[6])  # Ethanol exchange
    model3.reactions.get_by_id("y002056").bounds = (x[7], x[7])  # Succinate exchange
    model3.reactions.get_by_id("exchange_itaconate").bounds = (x[8],x[8])  # itaconate exchange
    solution1 = model3.optimize()
    print(x[0])
    print(solution1.objective_value)
    print(solution1.fluxes["y001663"])
    fluxes = solution1.fluxes.tolist()
    all_fluxes[x[0]] = fluxes

all_fluxes.to_excel("result/fluxes_for_all_condition.xlsx")






all_fluxes = pd.DataFrame({"rxnID": all_rxn})

for i, x in all_data.iterrows():
    #print(i)
    #print(i, x)
    qs = x[1]
    qco2 = x[2]
    qo2 = x[3]
    qcit = x[4]
    qpyr = x[5]
    qeth = x[6]
    qsucc = x[7]
    qita = x[8]
    growth = x[16]
    if x[14] == 'Ammonium Sulphate':
        model3.reactions.get_by_id("y001654").bounds = (-1000, 0)  # ammonia
        model3.reactions.get_by_id("y002091").bounds = (0, 0)  # urea
        model3.objective = {model3.reactions.y001714: 1}
        model3.reactions.get_by_id("y001714").bounds = (-1000, 0)  # D-glucose exchange
        model3.reactions.get_by_id("y002111").bounds = (growth, growth)  # growth exchange
        # model3.reactions.get_by_id("y001714").bounds = (qs*(-1), 0)  # D-glucose exchange
        model3.reactions.get_by_id("y001663").bounds = (qco2, qco2)  # bicarbonate exchange
        # model3.reactions.get_by_id("y001992").bounds = (qo2*(-1), qo2*(-1))  # oxygen exchange
        # by-product
        model3.reactions.get_by_id("y001793").bounds = (0, 0)  # formate exchange
        model3.reactions.get_by_id("y001634").bounds = (0, 0)  # acetate exchange
        model3.reactions.get_by_id("y001687").bounds = (x[4], x[4])  # citrate exchange
        model3.reactions.get_by_id("y002033").bounds = (x[5], x[5])  # Pyruvate exchange
        model3.reactions.get_by_id("y001761").bounds = (x[6], x[6])  # Ethanol exchange
        model3.reactions.get_by_id("y002056").bounds = (x[7], x[7])  # Succinate exchange
        model3.reactions.get_by_id("exchange_itaconate").bounds = (x[8], x[8])  # itaconate exchange
        solution1 = model3.optimize()
    else:
        model3.reactions.get_by_id("y001654").bounds = (0, 0) # ammonia
        model3.reactions.get_by_id("y002091").bounds = (-1000, 0) # urea
        solution1 = model3.optimize()

    fluxes = solution1.fluxes.tolist()
    all_fluxes[x[0]] = fluxes

all_fluxes.to_excel("result/fluxes_for_all_condition.xlsx")

