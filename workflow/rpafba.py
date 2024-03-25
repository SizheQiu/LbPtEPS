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
        out_medium[k] = 0.2 # maximum amino acid uptake rate
    for EX_m in (vitamins+DNA_materials):
        out_medium[EX_m] = 10
    for EX_m in others:
        out_medium[EX_m] = 1000
    return out_medium


def load_model_params(path):
    table = pd.read_csv(path)
    params = {'FpHmin':{} }
    A_table = table[table['Param']=='A'].reset_index().drop(['index'],axis=1)
    params['A'] = { A_table['Rxn'][i]: float(A_table['Value'][i]) for i in range(len(A_table.index)) }
    klach_table = table[table['Param']=='klach'].reset_index().drop(['index'],axis=1)
    params['klach'] = { klach_table['Rxn'][i]: float(klach_table['Value'][i]) for i in range(len(klach_table.index)) }
    Amin_table = table[table['Param']=='Amin'].reset_index().drop(['index'],axis=1)
    params['Amin'] = { Amin_table['Rxn'][i]: float(Amin_table['Value'][i]) for i in range(len(Amin_table.index)) }
    Fmin_table = table[table['Param']=='FpHmin'].reset_index().drop(['index'],axis=1)
    params['FpHmin'] = { Fmin_table['Rxn'][i]: float(Fmin_table['Value'][i]) for i in range(len(Fmin_table.index)) }
    return params


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



def get_FpH(rxn,  pH, ApH_table, FpHmin_dict):
    #update enzyme acticity based on pH.
    Fmin = FpHmin_dict[rxn]
    temp_pd = ApH_table[ApH_table['Enzyme']==rxn].reset_index().drop(['index'],axis=1)
    rxn, func_name, k1,k2 = list(ApH_table.iloc[0])
    if func_name == 'cubic':
        F = k1*(pH-k2)**3
    else:
        F = 1/(1+np.exp(-k1*(pH-k2)))
    return max(F, Fmin)
            

def set_LpPA( model, ptot, params, lac_con, pH, ApH_table ):
    sigma = 0.5
    #compute u%
    r0, r1, kpH, k1 = 0.0876, 0.07953, 44.9457 , 5.0344;
    u_ratio= r0+r1/( 1+np.exp(kpH*(pH-k1)) )
    
    A_dict = params['A']
    klachs = params['klach']
    Amin_dict = params['Amin']
    FpHmin_dict = params['FpHmin'] 
    #update LacH inhbition on carbon source uptake
    A_GLCpts = LacIn( A_dict['GLCpts'], Amin_dict['GLCpts'], lac_con, klachs['GLCpts'] )
    A_MANpts = LacIn( A_dict['MANpts'], Amin_dict['MANpts'], lac_con, klachs['MANpts'] )
    A_LCTSt = LacIn( A_dict['LCTSt'], Amin_dict['LCTSt'], lac_con, klachs['LCTSt'] )
    # A and T sectors
    p_sector = model.reactions.biomass_LPL60.flux_expression/A_dict['biomass_LPL60'] +\
           model.reactions.GLCpts.flux_expression/A_GLCpts  +\
           model.reactions.MANpts.flux_expression/A_MANpts +\
           model.reactions.LCTSt.flux_expression/A_LCTSt +\
           model.reactions.EX_ac_e.flux_expression/(sigma*A_dict['EX_ac_e']) + \
           model.reactions.EX_lac_L_e.flux_expression/(sigma*A_dict['EX_lac_L_e']) +\
           model.reactions.EX_lac_D_e.flux_expression/(sigma*A_dict['EX_lac_D_e'])
    # C sector 
    for k in a_dict.keys():
        if k not in ['EX_ac_e','EX_lac_L_e', 'EX_lac_D_e','biomass_LPL60', 'GLCpts','MANpts','LCTSt',\
                     'MANT_EPS', 'GLCT_EPS', 'GALT_EPS', 'WZX']:
            if k in FpHmin_dict.keys():
                F = get_FpH( k, pH, ApH_table, FpHmin_dict) # add pH inhibition on A
                p_sector  = p_sector  + model.reactions.get_by_id(k).flux_expression/( sigma*F*A_dict[k] )
            else:
                p_sector  = p_sector  + model.reactions.get_by_id(k).flux_expression/( sigma*A_dict[k] )
            
    # U sector
    F_GT = get_FpH( 'GT', pH, ApH_table, FpHmin_dict)
    u_sector= model.reactions.MANT_EPS.flux_expression/( sigma*F_GT*A_dict['MANT_EPS'] ) +\
            model.reactions.GLCT_EPS.flux_expression/( sigma*F_GT*A_dict['GLCT_EPS'] ) +\
            model.reactions.GALT_EPS.flux_expression/( sigma*F_GT*A_dict['GALT_EPS'] )+\
            model.reactions.WZX.flux_expression/( sigma*A_dict['WZX'] )
              
    PA = model.problem.Constraint( expression = p_sector + u_sector ,name = 'PA', lb= 0, ub = 0.5*ptot)
    pHrPA = model.problem.Constraint( expression = u_sector - u_ratio * ( p_sector + u_sector ),
                                                     name = 'pHrPA', lb= 0, ub = 0 )
    
    model.add_cons_vars([ PA, pHrPA ])
    
    return p_sector, u_sector



