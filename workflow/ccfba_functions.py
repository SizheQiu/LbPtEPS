import cobra
from cobra.io import read_sbml_model
from cobra import Model, Reaction, Metabolite
import numpy as np



cofactors = ['nad','atp','gtp','ctp','utp','ttp',
             'fad','adp','gdp','cdp','udp','tdp','coa']
macro = ['quino','cytoch','heme']


def checkMet( met ):
    """
    check whether metabolite contributes to carbon flux
    """
    for cof in cofactors:
        if cof in met.id:
            return False
    for ma in macro:
        if ma in met.name:
            return False
    return True

def addBound(model, rxn_id, num_atom, f_total ):
    """
    apply ub and lb for the reaction.
    """
    if num_atom == 0:
        return 0
    else:
        temp_rxn = model.reactions.query( rxn_id )
        if temp_rxn.reversibility == True:
            temp_rxn.lower_bound = - (f_total/num_atom)
            temp_rxn.upper_bound = f_total/num_atom
        else:
            temp_rxn.upper_bound = f_total/num_atom
            
    return 1


def getCNCount( rxn ):
    """
    count C,N atoms flowing through the reaction.
    """
    C,N = 0,0
    for met in rxn.reactants:
        if checkMet( met ):
            coeff =  abs(rxn.get_coefficient(met.id))
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
        
              

def ccfba(model):
    
    C_in, N_in = getCN_influx(model)['C'], getCN_influx(model)['N']
    
    for rxn in model.reactions:
        if 'EX_' in rxn.id or 'sink_' in rxn.id or 'Growth' in rxn.id:
            continue
        else:
            CN_count = getCNCount(rxn)
            addBound(model, rxn.id, CN_count['C'], C_in )
            addBound(model, rxn.id, CN_count['N'], N_in )
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            



