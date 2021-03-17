# -*- coding: utf-8 -*-
# -*- python 3 -*-
# -*- hongzhong Lu -*-


from cobra.io import read_sbml_model, load_matlab_model, write_sbml_model
import cobra
cobra_config = cobra.Configuration()
cobra_config.solver = "glpk_exact"
import os
import sys
import pprint
import pandas as pd

# model = read_sbml_model("../data/iYLI647_update2.xml") # can't read xml file now
model = load_matlab_model("model/iYLI647_update2.mat")

# nitrogen limitation
model.objective = 'R_r156'
model.reactions.get_by_id("R_EX_glc(e)").bounds = (-0.67,0)
model.reactions.get_by_id("R_EX_o2(e)").bounds = (-2,0)

solution = model.optimize()
print(model.summary())
fluxes = pd.DataFrame(solution.fluxes)
fluxes.to_csv("result/iYLI647.csv")

# test using the new data from Fujing
model = load_matlab_model("model/iYLI647_update2.mat")
# strain 1
model.reactions.get_by_id("R_EX_glc(e)").bounds = (-1.04, 0)
#model.reactions.get_by_id("R_EX_o2(e)").bounds = (-2,0) # current no data
model.reactions.get_by_id("R_EX_co2(e)").bounds = (2.06, 2.06)
solution2 = model.optimize()

# strain 2
model.reactions.get_by_id("R_EX_glc(e)").bounds = (-1.237, 0)
#model.reactions.get_by_id("R_EX_o2(e)").bounds = (-2,0) # current no data
model.reactions.get_by_id("R_EX_co2(e)").bounds = (3.00, 3.00)
solution2 = model.optimize()


# test using Ed data
model = load_matlab_model("model/iYLI647_update2.mat")
model.reactions.get_by_id("R_EX_glc(e)").bounds = (-0.61, 0)
#model.reactions.get_by_id("R_EX_o2(e)").bounds = (-2,0) # current no data
model.reactions.get_by_id("R_EX_co2(e)").bounds = (1.5, 1.5)
solution3 = model.optimize()

