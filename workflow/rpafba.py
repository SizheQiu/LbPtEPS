import cobra
from cobra.io import read_sbml_model
from cobra import Model, Reaction, Metabolite
import numpy as np
'''
Regulatory proteome constrained flux balance analysis of L.plantarum
'''


def init_MRSmedium():
    out_medium = {}
    vitamins = ['EX_btn_e','EX_pnto_R_e', 'EX_ribflv_e', 'EX_thm_e', 'EX_fol_e', 'EX_nac_e','EX_pydam_e']
    DNA_materials = ['EX_ade_e', 'EX_gua_e', 'EX_xan_e', 'EX_ura_e','EX_thymd_e']
    others = ['EX_mn2_e','EX_so4_e', 'EX_h2o_e', 'EX_h_e', 'EX_pi_e', 'EX_nh4_e']
    out_medium['EX_glc_e'] = 100
    for k in ['EX_his_L_e', 'EX_ile_L_e', 'EX_leu_L_e', 'EX_lys_L_e', 'EX_met_L_e', 
                     'EX_phe_L_e', 'EX_thr_L_e', 'EX_trp_L_e', 'EX_val_L_e', 'EX_ala_L_e', 
                     'EX_arg_L_e', 'EX_asn_L_e', 'EX_asp_L_e', 'EX_cys_L_e', 'EX_gln_L_e', 
                     'EX_glu_L_e', 'EX_gly_e', 'EX_pro_L_e', 'EX_ser_L_e', 'EX_tyr_L_e']:
        out_medium[k] = 0.2
    for EX_m in (vitamins+DNA_materials):
        out_medium[EX_m] = 10
    for EX_m in others:
        out_medium[EX_m] = 1000
    return out_medium

def approx_pH(lac_con, c1, c2):
    #approximate pH change in lactic acid fermentation
    return c1*lac_con + c2

def approx_pH(x):
    k1=4438.85;k2=7.62;
    #approximate pH change in MRS medium
    return np.log(k1/(x+k2) )

def LacIn( v_max,v_min, lac_con, klach ):
    pH = approx_pH(lac_con)
    lach = lac_con/(10**(pH-3.86))
    temp_value = abs(v_max)*exp( -klach*lach )
    if v_max > 0:
        out = max(temp_value, abs(v_min) )
    else:
        out = -1*max( temp_value, abs(v_min) )
    return out



def update_activity(A_dict, pH, lac_con, ApH_table):
    #update enzyme acticity based on pH.
    updated_As = {}
    lach = lac_con/(10**(pH-3.86))
    updated_As['GLCpts'] = abs(A_dict['GLCpts'])*exp( -klach*lach )
    for i in range(len(ApH_table.index)):
        rxn, func_name, k1,k2 = list(ApH_table.iloc[i])
        A0=A_dict[rxn]
        if func_name == 'cubic':
            updated_As[rxn] = A0 * ( k1*(pH-k2)**3 )
        else:
            updated_As[rxn] = A0 * ( 1/(1+np.exp(-k1*(pH-k2))) )
    return updated_As
            

def set_LpPA( model, ptot, params, profile, pH, a_dict, amin_dict ):
    sigma = 0.5
    klach = params['klach']
    #update LacH inhbition on carbon source uptake
    lac_con = profile['lac_D_e']
    A_GLCpts = LacIn( a_dict['GLCpts'], amin_dict['GLCpts'], lac_con, klach )
    A_MANpts = LacIn( a_dict['MANpts'],amin_dict['MANpts'], lac_con, klach )
    A_LCTSt = LacIn( a_dict['LCTSt'],amin_dict['LCTSt'], lac_con, klach )
    # A and T sectors
    expr = model.reactions.biomass_LPL60.flux_expression/a_dict['biomass_LPL60'] +\
           model.reactions.GLCpts.flux_expression/(-1*A_GLCpts)  +\
           model.reactions.MANpts.flux_expression/(-1*A_MANpts) +\
           model.reactions.LCTSt.flux_expression/(-1*A_LCTSt) +\
           model.reactions.EX_ac_e.flux_expression/(sigma*a_dict['EX_ac_e']) + \
           model.reactions.EX_lac_L_e.flux_expression/(sigma*a_dict['EX_lac_L_e']) +\
           model.reactions.EX_lac_D_e.flux_expression/(sigma*a_dict['EX_lac_D_e'])
    # C and U sectors     
    for k in a_dict.keys():
        if ( 'EX_' not in k ) and ( 'biomass' not in k ):
            expr = expr + model.reactions.get_by_id(k).flux_expression/( sigma*a_dict[k] )
              
    PA = model.problem.Constraint( expression = expr,name = 'PA', lb= 0, ub = 0.5*ptot)
    model.add_cons_vars([ PA ])
    
    return expr

def set_LppHPA( model, ptot, params, profile, pH, a_dict, amin_dict ):
    sigma=0.5
    r0, r1, kpH, k1 = 0.0876, 0.07953, 44.9457 , 5.0344;
    u_ratio= r0+r1/( 1+np.exp(kpH*(pH-k1)) )
    
    klach = params['klach']
    #update LacH inhbition on carbon source uptake
    lac_con = profile['lac_D_e']
    A_GLCpts = LacIn( a_dict['GLCpts'], amin_dict['GLCpts'], lac_con, klach )
    A_MANpts = LacIn( a_dict['MANpts'],amin_dict['MANpts'], lac_con, klach )
    A_LCTSt = LacIn( a_dict['LCTSt'],amin_dict['LCTSt'], lac_con, klach )
    
    p_sector=model.reactions.biomass_LPL60.flux_expression/a_dict['biomass_LPL60'] +\
           model.reactions.EX_glc_e.flux_expression/(-1*a_dict['EX_glc_e'])  +\
           model.reactions.EX_ac_e.flux_expression/(sigma*a_dict['EX_ac_e']) + \
           model.reactions.EX_lac_L_e.flux_expression/(sigma*a_dict['EX_lac_L_e'])
    
    for k in a_dict.keys():
        p_sector = p_sector + model.reactions.get_by_id(k).flux_expression/( sigma*a_dict[k] )
        
    u_sector= model.reactions.MANT_EPS.flux_expression/( sigma*a_dict['MANT_EPS'] ) +\
            model.reactions.GLCT_EPS.flux_expression/( sigma*a_dict['GLCT_EPS'] ) +\
            model.reactions.GALT_EPS.flux_expression/( sigma*a_dict['GALT_EPS'] )+\
            model.reactions.WZX.flux_expression/( sigma*a_dict['WZX'] )
               
    
    pHPA = model.problem.Constraint( expression = u_sector - u_ratio * ( p_sector + u_sector ),
                                                     name = 'pHPA', lb= 0, ub = 0 )
    model.add_cons_vars([ pHPA ])
    return u_sector - ratio * ( p_sector + u_sector )



