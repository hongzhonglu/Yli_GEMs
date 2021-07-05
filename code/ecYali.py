# -*- coding: utf-8 -*-
# -*- python 3 -*-
# -*- hongzhong Lu -*-

# Find the gene targets whose deletion could reduce lipid production
from cobra.io import read_sbml_model, load_matlab_model, write_sbml_model
import cobra
cobra_config = cobra.Configuration()
from cobra.manipulation import *
import pandas as pd
from cobra import Reaction
# eciYali can't predict the normal growth!!
# model0 = read_sbml_model("eciYali/model/eciYali.xml") # the reaction id should not contain white space!!


model0 = read_sbml_model("eciYali/model/eciYali_batch.xml") # the reaction id should not contain white space!!
solution1 = model0.optimize()

# function
def ecYliMinimalMedia(model, glucose_uptake=5):
    """
    This function is used to define a simple media for ecYeast
    :param model:
    :return: a model with the defined the minimal media
    """

    # test
    # model = model0

    rxnID = []
    rxnName = []
    for i, x in enumerate(model.reactions):
        rxnID.append(x.id)
        rxnName.append(x.name)

    exchange_rxn =[x for x, y in zip(rxnID, rxnName) if 'prot_' not in x and 'exchange' in y]
    uptake_rxn =[x for x in exchange_rxn if "_REV" in x]
    # first block any uptake
    for i, x in enumerate(uptake_rxn):
        rxn0 = uptake_rxn[i]
        print(rxn0)
        model.reactions.get_by_id(rxn0).upper_bound = 0

    # Block O2 and glucose production(avoids multiple solutions):
    model.reactions.get_by_id("y001992").upper_bound = 0  # 'oxygen exchange'
    model.reactions.get_by_id("y001714").upper_bound = 0  # 'D-glucose exchange'
    model.reactions.get_by_id("y001714_REV").upper_bound = glucose_uptake  # 'D-glucose exchange'

    #Allow uptake of essential components
    model.reactions.get_by_id("y001654_REV").upper_bound = 10000 #'ammonium exchange'
    model.reactions.get_by_id("y002100_REV").upper_bound = 10000 #'water exchange'
    model.reactions.get_by_id("y001861_REV").upper_bound = 10000 #'iron(2+) exchange'
    model.reactions.get_by_id("y001992_REV").upper_bound = 10000 #'oxygen exchange'
    model.reactions.get_by_id("y002005_REV").upper_bound = 10000 #'phosphate exchange'
    model.reactions.get_by_id("y002060_REV").upper_bound = 10000 #'sulphate exchange'
    model.reactions.get_by_id("y001832_REV").upper_bound = 10000 #'H+ exchange'
    model.reactions.get_by_id("y001671_REV").upper_bound = 10000 #Biotin . exchange
    model.reactions.get_by_id("y001548_REV").upper_bound = 10000 #pantothenate exchange
    model.reactions.get_by_id("y001967_REV").upper_bound = 10000 #Niconitate exchange
    model.reactions.get_by_id("y001947_REV").upper_bound = 10000 #Myo-inositol
    model.reactions.get_by_id("y002067_REV").upper_bound = 10000 #Thiamin (1+) exchange
    model.reactions.get_by_id("y002028_REV").upper_bound = 10000 #Pyridoxine exchange
    model.reactions.get_by_id("y001604_REV").upper_bound = 10000 #Aminobenzoic acid
    # Block bicarbonate uptake
    model.reactions.get_by_id("y001663").upper_bound = 0  # 'bicarbonate uptake'

    # solution1 = model0.optimize()

    # check the result
    # result_df = pd.DataFrame({"rxnid": rxnID, "rxnName": rxnName, "fluxes": solution1.fluxes})
    # result_df.to_excel("result/flux_ecModel.xlsx")

    return model

model1 = ecYliMinimalMedia(model=model0)
solution2 = model1.optimize()



def getModelWithRemoveGene (model0, gene_remove0):
    model1 = model0.copy()
    remove_genes(model1, gene_remove0, remove_reactions=True) # this script can't be used
    return model1

def solveEcModel(model, obj_id, target_id, glucose_uptake ='y001714_REV'):
    with model:
        model.objective = obj_id  # we put it in the mutant model
        solution = model.optimize()  # we put it in the mutant model
        # mutModel.summary() # we put it in the mutant model
        solution_value0 = solution.objective_value  # we put it in the mutant model
        model.reactions.get_by_id(obj_id).bounds = (0.999 * solution_value0, solution_value0)  # we put it in the mutant model
        model.objective = {model.reactions.prot_pool_exchange: -1}
        solution2 = model.optimize()
        yield_p = solution2.fluxes[target_id]/solution2.fluxes[glucose_uptake]
    return yield_p


# prepare gene ID
geneID = []
for x in model0.genes:
    print(x.id)
    geneID.append(x.id)

# using pfba
yield_p_ref0 =[]
yield_p_mutant0 =[]
for i, x in enumerate(geneID):
    # test
    print(i)
    gene_remove = []
    gene_remove.append(x)
    mutModel = getModelWithRemoveGene(model0=model1, gene_remove0=gene_remove)

    # solve the reference strain
    try:
        yield_p_ref = solveEcModel(model=model1, obj_id='xLIPID', target_id='xLIPID')
        print(yield_p_ref)
        yield_p_ref0.append(yield_p_ref)
    except:
        yield_p_ref = None
        yield_p_ref0.append(yield_p_ref)

    # solve the mutant strain
    try:
        yield_p_mutant = solveEcModel(model=mutModel, obj_id='xLIPID', target_id='xLIPID')
        print(yield_p_mutant)
        yield_p_mutant0.append(yield_p_mutant)
    except:
        yield_p_mutant = None
        yield_p_mutant0.append(yield_p_mutant)





# evalulate the effect of gene deletion on the growth

# reaction 1
dict1 = {model1.metabolites.get_by_id('s_0434[c]'): -1,
         model1.metabolites.get_by_id('m1726[c]'):-1,
         model1.metabolites.get_by_id('carbohydrate[c]'):-1,
         model1.metabolites.get_by_id('RNA[c]'): -1,
         model1.metabolites.get_by_id('DNA[c]'): -1,
         model1.metabolites.get_by_id('s_0394[c]'):1,
         model1.metabolites.get_by_id('s_1322[c]'):1,
         model1.metabolites.get_by_id('s_0450[c]'):1
        }

# using the general procedures to add the new function
reaction = Reaction('new1')
reaction.name = 'xBIOMASS_no_lipid'
reaction.subsystem = 'new_added'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.EC = ''
reaction.add_metabolites(dict1)
reaction.gene_reaction_rule = ''
model1.add_reactions([reaction])



# using pfba
yield_x_ref0 =[]
yield_x_mutant0 =[]
for i, x in enumerate(geneID):
    # test
    print(i)
    gene_remove = []
    gene_remove.append(x)
    mutModel = getModelWithRemoveGene(model0=model1, gene_remove0=gene_remove)

    # solve the reference strain
    try:
        yield_x_ref = solveEcModel(model=model1, obj_id='new1', target_id='new1')
        print(yield_x_ref)
        yield_x_ref0.append(yield_x_ref)
    except:
        yield_x_ref = None
        yield_x_ref0.append(yield_x_ref)

    # solve the mutant strain
    try:
        yield_x_mutant = solveEcModel(model=mutModel, obj_id='new1', target_id='new1')
        print(yield_x_mutant)
        yield_x_mutant0.append(yield_x_mutant)
    except:
        yield_x_mutant = None
        yield_x_mutant0.append(yield_x_mutant)


pfba_result = pd.DataFrame({'gene':geneID, 'yield_p_ref0': yield_p_ref0, 'yield_p_mutant0': yield_p_mutant0, 'yield_x_ref0': yield_x_ref0, 'yield_x_mutant0': yield_x_mutant0})
pfba_result.to_excel("result/gene_deletion_prediction.xlsx")



