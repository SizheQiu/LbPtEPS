import cobra
from cobra.io import read_sbml_model
from cobra import Model, Reaction, Metabolite
import numpy as np



cofactors = ['nad','atp','gtp','ctp','utp','ttp',
             'fad','adp','gdp','cdp','udp','tdp']
macro = ['quino','cytoch','heme']
non_met = [ 'EX_','sink_', 'DM_', 'Growth','biomass'  ]

def isCofactor( met ):
    for cof in cofactors:
        if cof in met.id:
            return True
    for ma in macro:
        if ma in met.name:
            return True
    return False

def isCoA( met ):
    if 'coa' in met.id:
        return True
    return False
    
def isACP( met ):
    if 'ACP' in met.id:
        return True
    return False

def isNonMet( rxn ):
    for prefix in non_met:
        if prefix in rxn.id:
            return True
    return False
    

def getCNCount( rxn ):
    """
    count C,N atoms flowing through the reaction.
    """
    C,N = 0,0
    for met in rxn.reactants:
        if ! (isCofactor( met ) ):
            coeff =  abs(rxn.get_coefficient(met.id))
            if isCoA( met ):
                C += int(coeff)*( met.elements['C'] - 21 )
                N += int(coeff)*( met.elements['N'] - 7 )
            elif isACP( met ):
                C += int(coeff)*( met.elements['C'] - 11 )
                N += int(coeff)*( met.elements['N'] - 2 )
            else:
                C += int(coeff)*met.elements['C']
                N += int(coeff)*met.elements['N']
                
    return {'C':C,'N':N}
            

def getCN_influx(model):
    C_in,N_in = 0,0
    for exchange_id in model.medium:
        ex_rxn = model.reactions.get_by_id(exchange_id)
        CN_count = getCNCount(ex_rxn)
        C_in += CN_count['C']*model.medium[exchange_id]
        N_in += CN_count['N']*model.medium[exchange_id]
    return {'C':C_in,'N':N_in}
        

    
def addBound(model, rxn_id, num_atom, f_total ):
    """
    apply ub and lb for the reaction.
    """
    if num_atom == 0:
        return 0

    temp_rxn = model.reactions.query( rxn_id )
    bound = f_total/num_atom
    old_lb, old_ub = temp_rxn.lower_bound, temp_rxn.upper_bound
    if temp_rxn.reversibility == True:
        temp_rxn.lower_bound = max( - bound, old_lb )
        temp_rxn.upper_bound = min( bound, old_ub )
    else:
        temp_rxn.upper_bound = min( bound, old_ub )
            
    return 1
              

def ccfba(model):
    
    C_in, N_in = getCN_influx(model)['C'], getCN_influx(model)['N']
    for rxn in model.reactions:
        if  isNonMet(rxn):
            continue
        else:
            CN_count = getCNCount(rxn)
            addBound(model, rxn.id, CN_count['C'], C_in )
            addBound(model, rxn.id, CN_count['N'], N_in )
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            



