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

os.chdir('/Users/luho/PycharmProjects/model/cobrapy/code')
sys.path.append(r"/Users/luho/PycharmProjects/model/cobrapy/code")
pprint.pprint(sys.path)
# model = read_sbml_model("../data/iYLI647_update2.xml") # can't read xml file now
model = load_matlab_model("/Users/luho/PycharmProjects/model/cobrapy/data/iYLI647_update2.mat")


# nitrogen limitation
model.objective = 'R_r156'
model.reactions.get_by_id("R_EX_glc(e)").bounds = (-0.67,0)
model.reactions.get_by_id("R_EX_o2(e)").bounds = (-2,0)

solution = model.optimize()
print(model.summary())
fluxes = pd.DataFrame(solution.fluxes)
fluxes.to_csv("../result/iYLI647.csv")