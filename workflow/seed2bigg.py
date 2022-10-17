import cobra
from cobra.io import read_sbml_model,load_json_model
from cobra import Model, Reaction, Metabolite
import pandas as pd
import numpy as np
import json
import re
import pickle
from six import string_types


def unify_ex( model ):
    for met in model.metabolites.query('_b'):
        met_id = met.id.replace('_b', '')
        rxn_ids = [ rxn.id for rxn in met.reactions ]
        if 'EX_'+met_id+'_e0' not in rxn_ids:
            temp_met = model.metabolites.get_by_id(met_id+'_b')
            temp_met.id = met_id+'_e0'
            temp_met.name = str(temp_met.name).replace('_b', '').replace('_c0', '_e')
            temp_met.compartment = 'e0'
            temp_rxn = model.reactions.get_by_id('EX_'+met_id+'_b')
            temp_rxn.id = 'EX_'+met_id+'_e0'
            temp_rxn.name = 'EX_'+met_id+'_e0' 
        else:
            model.remove_metabolites([model.metabolites.get_by_id(met_id+'_b')])
            model.remove_reactions([model.reactions.get_by_id('EX_'+met_id+'_b')])
    return model
    

def update_met( model, path2map ):
    '''
    Update metabolites'IDs and names from seed id to bigg id.
    '''
    map_met = pd.read_csv(path2map,index_col='SEED')
    map_compt = {'e0':'e', 'c0':'c'}
    updates = {}
    for met in model.metabolites:
        id0 = met.id.replace('_'+met.compartment, '')
        if id0 in map_met.index:
            id1 = map_met.loc[id0,'nocompt']
            if not isinstance(id1,string_types):
                id1 = id1.values[0]
            updates[met.id] = id1 + '_' + map_compt[met.compartment]
            met.id = id1 + '_' + map_compt[met.compartment]
            met.name = list(map_met[map_met['nocompt']==id1]['name'])[0]
        if met.compartment == 'e0':
            met.compartment = 'C_e'
        if met.compartment == 'c0':
            met.compartment = 'C_c'
    model.compartments['C_c']='cytosol'
    model.compartments['C_e']='extracellular space'
    model.repair()
    return updates

def update_rxn( model, path2map, path2stoich ):
    '''
    Update reaction IDs and Names by macthing stoichiometry and the map between bigg and modelseed.
    '''
    map_rxn = pd.read_csv(path2map,index_col='SEED')
    with open( path2stoich , 'rb') as f:
        udb_rxn_dict = pickle.load(f)
    rxn_updates = {}
    for rxn in model.reactions:
        if 'EX_' in rxn.id:
            if len(rxn.products) < 1:
                rxn.id = 'EX_' + rxn.reactants[0].id
                continue
            
        rxn_id0 = str(rxn.id).replace('_c0', '')
        stoich_str = {k.id:v for k,v in rxn.metabolites.items()}
        rev_stoich_str = {k.id:-v for k,v in rxn.metabolites.items() }
        if stoich_str in udb_rxn_dict.values():# update reaction id by stoichiometry first
            rxn_id1 = [k for k,v in udb_rxn_dict.items() if stoich_str==v][0]
            rxn_updates[str(rxn.id)] = rxn_id1
            rxn.id = rxn_id1
        elif rev_stoich_str in udb_rxn_dict.values():
            rxn_id1 = [k for k,v in udb_rxn_dict.items() if rev_stoich_str==v][0]
            rxn_updates[str(rxn.id)] = rxn_id1
            rxn.id = rxn_id1  
#         elif rxn_id0 in map_rxn.index:
#             rxn_id1 = map_rxn.loc[rxn_id0,'bigg']
#             if not isinstance(rxn_id1,string_types):
#                 rxn_id1 = rxn_id1.values[0]
#             rxn_updates[str(rxn.id)] = rxn_id1
#             rxn.id = rxn_id1
            
    model.repair()
    return rxn_updates

def update_ex( model ):
    for rxn_id in model.medium.keys():
        temp_rxn = model.reactions.get_by_id(rxn_id)
        reactant_id = temp_rxn.reactants[0].id
        if rxn_id != 'EX_' + reactant_id:
            temp_rxn.id = 'EX_' + reactant_id
            temp_name.id = 'EX_' + reactant_id
    return model
            
    








