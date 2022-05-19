# -*- coding: utf-8 -*-
# -*- python 3 -*-
# -*- hongzhong Lu -*-

# Find the gene targets whose deletion could reduce lipid production
from cobra.io import read_sbml_model, load_matlab_model, write_sbml_model, load_json_model
import cobra
cobra_config = cobra.Configuration()
from cobra.manipulation import *
import pandas as pd
from cobra import Reaction
import numpy as np

import os
import sys

from cobra import Reaction
from cobra.flux_analysis import pfba


from cameo.strain_design.deterministic.flux_variability_based import FSEOF

model0 = load_matlab_model("/Users/xluhon/Documents/GitHub/GECKO2_simulations/ecModels/ecKmarx/ecKmarx_batch.mat") # the reaction id should not contain white space!!
solution1 = model0.optimize()

GEM_model = read_sbml_model("/Users/xluhon/Documents/GitHub/Kluyveromyces_marxianus-GEM/modelFiles/xml/Kluyveromyces_marxianus-GEM.xml") # the reaction id should not contain white space!!
solution2 = GEM_model.optimize()


biomass_rxn_id0 = 'r_1913' # biomass
test_rxn_id0 = 'r_1889'# succinate


def sorted_fseof(model, biomass_rxn_id, test_rxn_id):
    # Revert the model to its original state:
    model.reactions.get_by_id(biomass_rxn_id).lower_bound = 0
    model.reactions.get_by_id(test_rxn_id).lower_bound = 0
    model.objective = biomass_rxn_id
    # Run analysis
    fseof = FSEOF(model)
    fseof_result = fseof.run(target=model.reactions.get_by_id(test_rxn_id))
    fseof_df = fseof_result.data_frame
    # For each row, create a linear model with the test exchange as prediction, and store the slope of said model:
    fseof_df["slope"] = np.nan
    fseof_df["r2"] = np.nan
    for index, row in fseof_df.iterrows():
        if sum(row) == 0:
            fseof_df.loc[index,"slope"] = 0
        else:
            x = row.iloc[:-2]
            y = fseof_df.loc[test_rxn_id].iloc[:-2]
            A = np.vstack([x, np.ones(len(x))]).T
            m, c = np.linalg.lstsq(A, y, rcond=None)[0]
            resid = np.linalg.lstsq(A, y, rcond=None)[1]
            r2 = 1 - resid / (y.size * y.var())
            fseof_df.loc[index,"slope"] = m
            try:
                fseof_df.loc[index,"r2"] = r2
            except:
                fseof_df.loc[index,"r2"] = 0
    # Sort the dataframe by slope and print only rows with R2 > 0.5:
    fseof_df = fseof_df.sort_values(by=["slope"], ascending=False)
    print(fseof_df.loc[fseof_df.index != "test"].loc[fseof_df["r2"] > 0.5].iloc[:20, :])
    return fseof_df

# Run the FSEOF analysis for both models:
fseof_df = sorted_fseof(model=GEM_model, biomass_rxn_id=biomass_rxn_id0, test_rxn_id=test_rxn_id0)

ec_fseof_df = sorted_fseof(model0, biomass_rxn_id0, test_rxn_id0)

fseof_df.to_excel("result/fseof_df_result_for_Kmarx.xlsx")
ec_fseof_df.to_excel("result/ec_fseof_df_result_for_Kmarx.xlsx")


# print the rxn ID and gene association
rxnID =[]
gpr=[]
for rxn in model0.reactions:
    print(rxn.id, rxn.gene_reaction_rule)
    rxnID.append(rxn.id)
    gpr.append(rxn.gene_reaction_rule)

ID_mapping = pd.DataFrame({"rxnID":rxnID,"gene":gpr})
ID_mapping.to_excel("result/ID_mapping_for_Kmarx.xlsx")


###




