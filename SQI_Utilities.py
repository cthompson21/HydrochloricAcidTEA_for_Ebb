# -*- coding: utf-8 -*-
"""
Created on Sat May 25 09:23:48 2024

@author: cthompson
"""


import numpy as np
from matplotlib.pyplot import *
import pdb
from numpy.linalg import inv
import pandas as pd

from matplotlib.pyplot import *
from matplotlib.pyplot import *


df_columns =['Feed Flow', 'Feed TDS', 'Permeate Flow', 'Permeate TDS',
       'Concentrate Flow', 'Concentrate TDS', 
       
       
       'Operating Pressure', 
       'HP Pump Efficiency', 'Circ Pump Efficiency', 'No. of CIP Blocks',
       
       'Vessels per Rack', 'Area per Rack', 
       'Unit Cost of Modules', 'Module Lifetime', 'Unit Cost of Vessels',
       'Flux', '# of Modules', '# of Vessels', 
        'HP Pump Power', 'Circ Flow Rate (GPM)',
       'Circulation Pump Power', 'Pump Operating Cost',
       'Annual Cost of Modules', 'Cost of Module Set', 'Cost of Vessels',
       'High Pressure Pumps', 'Circulation Pumps', 'Piping',
       'Valves and Instruments', 'Heat Exchangers', 'Selectivity',
       'MaxOsmPressure_psi', 'Circ Flow Rate', 'Circ Pump Power', 'NumVessels',
       'Peripheral Hardware', 'Total System Cost',
       'Construction and Assembly']

def calc_maximum_selectivity(max_tmp_Pa, xf,r, Vs, i, T ):
    
    """calculates the maximum selectivity allowed for a membrane, given a maximum transmembrane pressure allowed,
    according to the Van't Hoff. equation
    
    max_tmp_Pa: maximum allowable transmembrane pressure in pascals
    xf: the mole fraction of solids in feed
    r: is the single stage recovry, calculated as flow_permeate/flow_feed
    
    Vs: the molar volume of the pure solvent
    i: the van't hoff parameters (number of ions a solute would split up into)
    
    """
    
    R=8.314
    
    
    phi=np.exp(max_tmp_Pa*Vs/(-i*R*T))
    
    

    max_selectivity = (1/(1/xf-1))*((phi-r*phi+r)/(phi-r*phi-1+r+xf) -1)
    if isinstance(max_selectivity, (np.ndarray,)):
        max_selectivity[max_selectivity<=0]=np.inf
    else:
        if max_selectivity<=0:
            max_selectivity=np.inf
    return max_selectivity


# def calc_recovery_to_meet_membrane_selectivity(max_tmp_Pa, xf,membrane_selectivity, Vs, i, T, permeate_pressure=0 ):
    
#     """calculates the maximum recovery allowed for a membrane, given a maximum transmembrane pressure allowed,
#     according to the Van't Hoff. equation
    
#     max_tmp_Pa: maximum allowable transmembrane pressure in pascals
#     xf: the mole fraction of solids in feed
#     r: is the single stage recovry, calculated as flow_permeate/flow_feed
    
#     Vs: the molar volume of the pure solvent
#     i: the van't hoff parameters (number of ions a solute would split up into)
    
#     """
#     """
#     z= max_selectivity/(1/(1/xf-1)) 
#     (phi-r*phi+r)/(phi-r*phi-1+r+xf)    = z+1
   
#     (phi-r*phi+r)   = (z+1)*(phi-r*phi-1+r+xf) 
#     phi + (1-phi)*r   = (z+1) * (phi + xf-1 +r*(1-phi))
    
#     phi -  (z+1) * (phi + xf-1) = (z+1)*(1-phi)*r - (1-phi)*r  =z*(1-phi)*r 
    
#     r=( phi -  (z+1) * (phi + xf-1))/(z*(1-phi))
    
#     """
    
    
#     R=8.314
#     phi=np.exp((max_tmp_Pa+permeate_pressure)*Vs/(-i*R*T))
#     z= membrane_selectivity/(1/(1/xf-1)) 
    
#     r=( phi -  (z+1) * (phi + xf-1))/(z*(1-phi))

#     return r

def calc_recovery(feed_tds, membrane_selectivity, MW_solute, MW_solvent, osmotic_pressure_max_Pa,
                  i=1, Vs = 5.556e-5):
    """ performs numeric calculatuion of the maximum allowable recovery for a given allowable osmotic pressure"""
    
    x_feed = feed_tds
    mass_feed = 1
    mass_feed_solute= feed_tds*1
    mass_feed_solvent = (1-feed_tds)*1
    mass_perm_solvent = 0
    mass_perm_solute = 0
    
    
    permeate_ratio =np.array([])
    permeate_conc = np.array([])
    feed_conc = np.array([])
    mass_left = np.array([])
    mass_right = np.array([])
    for repeat in range(50):
        """  """
        d_m = 0.01
        Z = membrane_selectivity*(1-x_feed)/x_feed
        x_perm  = 1/(Z+1)
        d_m_solute = x_perm*d_m
        d_m_solvent = (1-x_perm)*d_m
        mass_feed_solute -= d_m_solute
        mass_feed_solvent -= d_m_solvent
        mass_perm_solvent += d_m_solvent
        mass_perm_solute += d_m_solute
        
        x_feed = mass_feed_solute/(mass_feed_solute+mass_feed_solvent)
        permeate_ratio = np.append(permeate_ratio, mass_perm_solvent+mass_perm_solute)
        permeate_conc = np.append(permeate_conc, mass_perm_solute/(mass_perm_solute + mass_perm_solvent))
        feed_conc = np.append(feed_conc, x_feed)
        mass_right = np.append(mass_right, mass_perm_solute+mass_perm_solvent)
        mass_left = np.append(mass_left, mass_feed_solute + mass_feed_solvent)
        osm = calc_osmotic_pressure(feed_conc[-1], permeate_conc[-1], MW_solute, MW_solvent, i=i, Vs= Vs)
        if osm>osmotic_pressure_max_Pa:
            break
     
    
        
    permeate_tds = permeate_conc[-1]
    conc_tds = feed_conc[-1]
    permeate_ratio=permeate_ratio[-1]
   
  
    
    
    # figure()
    # plot(permeate_ratio, permeate_conc, label='permeate')
    
  
    # plot(permeate_ratio,feed_conc, label = 'feed')
    # plot(permeate_ratio,mass_left, label = 'left')
    # plot(permeate_ratio,mass_right, label = 'right')
    
    return permeate_tds, conc_tds, permeate_ratio, 
    







def solve_stage_feed_mixed_with_permeate(max_tmp_Pa,feed_flow, feed_tds ,
                                                    membrane_selectivity,
                                                    Vs,
                                                    i,
                                                    T,
                                                    MW_solute,
                                                    MW_solvent,
                                                    feed_bypass_portion, 
                                                    optimize_feed_bypass=True
                                                    
                                                    ):
    
   
    """
    solves a stage for a bypass flow -equipped stage, where some of the feed bypasses the membrane and is mixed with the permeate.
    
    max_tmp_Pa is the maximum absolute pressure allowed by the module in pascals
    feed_flow is the feed flow rate in any unit (usually GPM)
    feed_tds  is the mass fraction of feed that is solute
    membrane selectivity is the membrane selectiviey (should be greater than 1)
    Vs is the molar volume of the solvent
    i is the van't hoff coefficient
    MW_solute and solvnet are molecular weights in g/mol
    feed_bypass_portion is a number between 0 and 1 that gives the fractio nof the feed that is fed around to the permeate.
    
    
    outputs a tuple with various stage parameters
    
    the output feed flow and tds are different from the input feed flow and tds.  The output feed flow is the feed flow that reaches the membrane (does not bypass)
    the true permeate is the flow rate through the membrane, while the permeate is hte combined bypass flow and true permeate.
    concentrate is as usual
    
    """
    
    if np.any(np.isnan(feed_flow)):
        print('feed flow nan')
        pdb.set_trace()
    if np.any(np.isnan(feed_tds)):
        print('feed tds nan')
        pdb.set_trace()
    if isinstance(feed_tds, (list, np.ndarray)):
        if len(feed_tds)==1:
            feed_tds= feed_tds[0]
        
        
     #### first check for simple solution where oxmostic pressure is not limiting.  If it is possible to get above 10% recovery without a bypass, then we will accept that.  
     #### if not, we will go on to use the bypass solutoin
    
    permeate_tds, concentrate_tds, permeate_ratio,  = calc_recovery(feed_tds, membrane_selectivity, MW_solute, MW_solvent, max_tmp_Pa)
    if permeate_ratio>0.1:
        true_permeate_flow = permeate_ratio*feed_flow
        true_permeate_tds = permeate_tds
        permeate_flow = permeate_ratio*feed_flow
        concentrate_flow = feed_flow-true_permeate_flow
        bypass_ratio=0
        print('dont need to do bypass')
        return feed_flow, feed_tds, true_permeate_flow, true_permeate_tds, permeate_flow, permeate_tds, concentrate_flow, concentrate_tds
        
      
        
    
    last_solvent_removed=0.000001
    delta=1
    feed_flow_orig = float(feed_flow)
    cycle=0
    last_concentrate_tds = 0.001
    record = []
    while True: 
        cycle +=1
        
        permeate_ratio =np.arange(0.001,1,0.001)
        
        ### bypass here refers to the flow bypassing membrane
        bypass_flow = feed_bypass_portion*feed_flow_orig
        bypass_tds = feed_tds
        
        feed_flow = (1-feed_bypass_portion)*feed_flow_orig  ## feed_flow to membrane
        
        Z= (1-feed_tds)/feed_tds*membrane_selectivity
        true_permeate_tds=1/(Z+1)
        true_permeate_flow = feed_flow*permeate_ratio
        
        concentrate_flow=feed_flow-true_permeate_flow
        concentrate_tds= (feed_flow*feed_tds-true_permeate_flow*true_permeate_tds)/concentrate_flow

        permeate_tds = (bypass_flow*bypass_tds + true_permeate_flow* true_permeate_tds)/(bypass_flow + true_permeate_flow)
        
        osmotic_pressure = -calc_osmotic_pressure(concentrate_tds, permeate_tds,
                                  MW_solute, MW_solvent, i = i,Vs=Vs, T=T )
        
        to_del = np.append(np.where(np.isnan(osmotic_pressure))[0],np.where(osmotic_pressure>0)[0])
        osmotic_pressure = np.delete(osmotic_pressure, to_del)
        permeate_ratio = np.delete(permeate_ratio, to_del)
        
        if np.all(osmotic_pressure>max_tmp_Pa):
            print('osmotic pressure issure')
            pdb.set_trace() 
            return 0.001
        
        x=np.argmin(np.abs(osmotic_pressure+ max_tmp_Pa))

        permeate_ratio = permeate_ratio[x]

        true_permeate_flow = feed_flow*permeate_ratio
        Z= (1-feed_tds)/feed_tds*membrane_selectivity
        true_permeate_tds=1/(Z+1)
        
        concentrate_flow=feed_flow-true_permeate_flow
        concentrate_tds= (feed_flow*feed_tds-true_permeate_flow*true_permeate_tds)/concentrate_flow
        
        permeate_flow = true_permeate_flow+bypass_flow
        permeate_tds = (bypass_flow*bypass_tds + true_permeate_flow* true_permeate_tds)/(bypass_flow + true_permeate_flow)
        
        solvent_removed = true_permeate_flow*(1-true_permeate_tds)

        delta = (solvent_removed - last_solvent_removed)/last_solvent_removed
        last_solvent_removed = solvent_removed
        
        
        delta = (concentrate_tds-last_concentrate_tds)/last_concentrate_tds
        last_concentrate_tds = concentrate_tds*1
        
        record.append(delta)
        # print(feed_bypass_portion, delta)
        if cycle>25:
            
            print('did not find optimum bypass in 25 cycles. Last delta=', np.round(record[-2:],3), last_solvent_removed)
            break
        
        if np.abs(delta)<0.05 or solvent_removed<1e-6:
            break
        elif delta>0:
            feed_bypass_portion+=0.01
        elif delta<0:
            feed_bypass_portion-=0.01
        
    return feed_flow, feed_tds, true_permeate_flow, true_permeate_tds, permeate_flow, permeate_tds, concentrate_flow, concentrate_tds
    
    


def calc_osmotic_pressure(feed_conc_wtpercent, permeate_conc_wtpercent,
                          MW_solute, MW_solvent, i = 1,Vs=5.556e-05, T=298,
                          Vs_conc = None):
    
    """The feed concentration refers to the TDS on the high conc side.  It might actually correspond to concentrate TDS
    MW in g/mol
    Vs in m^3/mol
    
    i is the van't hoff coefficient
    T is temp in Kelvin
    """
    
    R=8.314
    
    if Vs_conc == None:
        Vs_conc = Vs
    
    x_solvent_permeate = (1-permeate_conc_wtpercent)/MW_solvent/((1-permeate_conc_wtpercent)/MW_solvent+i*permeate_conc_wtpercent/MW_solute)
    x_solvent_concentrate = (1-feed_conc_wtpercent)/MW_solvent/((1-feed_conc_wtpercent)/MW_solvent+i*feed_conc_wtpercent/MW_solute)

    osmotic_pressure = -i*R*T*np.log(x_solvent_concentrate/x_solvent_permeate)/Vs
    # osmotic_pressure = -i*R*T*np.log(x_solvent_concentrate)/Vs_conc+i*R*T*np.log(x_solvent_permeate)/Vs
    return osmotic_pressure



    



def solve_flows(df, feed_flow, feed_tds, feed_input_stage=None):
    
    """solves flows for a passback configuration"""
    
    N = len(df.index)
    p_over_f = df['Permeate Flow']/df['Feed Flow']
    F=np.zeros((N,))
    # pdb.set_trace()
    if feed_input_stage==None:
        # feed_input_stage = np.where((df['Feed TDS'].iloc[1:]>feed_tds).to_numpy() & (df['Feed TDS'].iloc[:-1]<feed_tds).to_numpy())[0]
        feed_input_stage=np.intersect1d(np.where(df['Feed TDS'].iloc[1:]>=feed_tds)[0],
                                    np.where(df['Feed TDS'].iloc[:]<feed_tds)[0])
        if len(feed_input_stage)==0:
            feed_input_stage=0
        else:
            feed_input_stage=feed_input_stage[0]
    elif isinstance(feed_input_stage, (int,)):
        pass
  
    F[feed_input_stage] = -feed_flow                        
    x=np.zeros((N,N))

    for i in np.arange(N):
        x[i,i]=-1/p_over_f[i]
        if i>0:
            x[i,i-1] = (1-p_over_f[i-1])/p_over_f[i-1]
    
    permeate_connections =df['Permeate Connections']
    for i in range(1,len(df.index)):
        permeate_cnct_stage = permeate_connections[i]
        x[permeate_cnct_stage, i]=1
        
    p = np.matmul(inv(x),F)
    f = p/p_over_f
    c=p*(1-p_over_f)/p_over_f
    

    return f,p,c,permeate_connections, feed_input_stage


def solve_conc(df, feed_flow, feed_tds,feed_input_stage=None, pause =False):
    
    
    """use same matrix eqations to calculate the solute flows in the system"""
    
    
    N = len(df.index)
    p_over_f = df['Permeate Flow']*df['Permeate TDS']/(df['Feed Flow']*df['Feed TDS'])
    F=np.zeros((N,))
    if pause:
        pdb.set_trace()
    
    if feed_input_stage==None:
        feed_input_stage=np.intersect1d(np.where(df['Feed TDS'].iloc[1:]>=feed_tds)[0],
                                        np.where(df['Feed TDS'].iloc[:]<feed_tds)[0])
        if len(feed_input_stage)==0:
            feed_input_stage=0
        else:
            feed_input_stage=feed_input_stage[0]
            
    elif isinstance(feed_input_stage, (int,)):
        pass
    
    # print('Feed onput', feed_input_stage)
    F[feed_input_stage] = -feed_flow*feed_tds                     
    x=np.zeros((N,N))

    for i in np.arange(N):
        x[i,i]=-1/p_over_f[i]  ### this is the permeat flow forward
        if i>0:
            x[i,i-1] = (1-p_over_f[i-1])/p_over_f[i-1]  ### this is hte concentrate from previous stage being passed in
    permeate_connections =df['Permeate Connections']
    for i in range(1,len(df.index)):
        permeate_cnct_stage = permeate_connections[i]          #### connects the higher permeates  
        x[permeate_cnct_stage, i]=1
        
    
    p = np.matmul(inv(x),F)
    
    f = p/p_over_f
    c=p*(1-p_over_f)/p_over_f
    
    # if np.any(p>df['Permeate Flow']):
    #     pdb.set_trace()
    # if np.any(f>df['Feed Flow']):
    #     pdb.set_trace()
    # if np.any(c>df['Concentrate Flow']):
    #     pdb.set_trace()
    
    ### convert to concentrations
    p=p/df['Permeate Flow']
    c=c/df['Concentrate Flow']
    f=f/df['Feed Flow']
   
    

    return f,p,c

if __name__ =="__main__":
    
   
    
    # mw= 265
    # x= calc_osmotic_pressure(0.4, 0.01, mw, 18, Vs=0.018/1000, T=373)/101325*14.7
    # print(x)
    
    # mw=20
    # x= calc_osmotic_pressure(0.5, 0.00, mw, 18, Vs=0.018/1000, i=2,T=373)/101325*14.7
    
    
    # print(x)
    
    
    # test = pd.DataFrame(columns = ['Feed Flow','Feed TDS','Concentrate Flow','Concentrate TDS', 'Permeate Flow','Permeate TDS'],
    #                     data=np.array([[200,0.25, 10,0.7,190, (200*0.25-10*0.7)/190 ],
    #                                   [400,0.1, 20,0.25,380, (200*0.1-10*0.25)/190 ],
    #                                   [100,0.02, 5,0.1,95, (200*0.02-10*0.1)/190 ],
    #                                   [50,0.004, 2.5,0.02,47.5, (200*0.004-10*0.02)/190 ],]))
    
    calc_recovery(0.019, 2, 108, 18 , 101325*100)
    