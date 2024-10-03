# -*- coding: utf-8 -*-
"""
Created on Thu Jan 18 11:48:44 2024

@author: cthompson


Still to do:  make the stage starts the same as the permeate connection values.  This should be the case otherwise, an extra pump is needed.

"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import *
import pdb
import os
import itertools
import seaborn as sns
import warnings
import copy

warnings.simplefilter(action='ignore', category=FutureWarning)
warnings.simplefilter(action='ignore', category=RuntimeWarning)
from distillation.solvers import solve_diagonal

from scipy.sparse import linalg, diags

pd.options.display.max_rows=100

pd.options.display.max_columns=20
pd.options.display.float_format = '{:,.5f}'.format
sns.set_palette("Paired")
sns.set(font_scale=1.5)
from sklearn.preprocessing import StandardScaler, MinMaxScaler
from sklearn.decomposition import PCA
save_file ='GreenMethanol_SDI_new.xlsx'
template_file = 'H:/Shared drives/Via Separations Main Folder/New Applications/108_sulfuric acid/0_Sulfuric acid TEA- Chris Mod.xlsx'
from SQI_Utilities import *

class Via(object):
    
    def __init__(self, feed_rate_mt_year, 
                 input_tds,
                 solute_recovery_target,
                 solvent_recovery_target,
                 concentrate_tds_target,
                 permeate_tds_target,

                 selectivity,  template=[], num_stages=[1,], flux=1, material='SS',
                 autooptimize = True,
                 levelized_cost_basis='concentrate', 
                 max_tmp_Pa = 1500/14.7*101325, ### 1000 psi in Pa
                 rho_solvent_kg_m3 =1000,
                 rho_solute_kg_m3 =1000,
                 MW_solute_g_mol = 32,
                 MW_solvent_g_mol = 18,
                 vantHoffi=1,
                 R=8.314,
                 T=298,
                 ):
        """ import feed rate in MT/year,
        input_tds as float between 0 and 1
        output_tds as float between 0 and 1
        selectivity as float greater than 1
        recovery  as float between 0 and 1
        
        
        In this version, the input tds, output tds target, membrane selectivity and recovery targets are used as inputs.  The optimize_self function finds a system design meeting the targets
        at lowest cost.  
        
        NOTE:  The choice of system_type is an erroneous feature.  In all cases, the permeate from a first pass should become the feed for the next pass.
        
        permeate recovery target is == permeate volume*(1-permeate_tds)/(feed_volume*(1-feed_tds))
        solute recovery target is == solute volume*solute_tds)/(feed_volume*feed_tds))
        """
        
        
        self.levelized_cost_basis=levelized_cost_basis
        self.material=material
        self.limit_osmotic_pressure=True
        
        self.selectivity = selectivity
        self.concentrate_tds_target = concentrate_tds_target     ### this is the target for output_tds
        self.permeate_tds_target = permeate_tds_target
        self.solute_recovery_target = solute_recovery_target
        self.solvent_recovery_target = solvent_recovery_target  ### this is the target for total recovery
        
        
        
        
        self.rho_solvent = rho_solvent_kg_m3
        self.rho_solute = rho_solute_kg_m3
        Vs_solvent = 0.098/1800
        Vs_solute = 0.098/1800
        Vs = (Vs_solvent + Vs_solute)/2
        
        self.max_tmp_Pa =max_tmp_Pa
        self.MW_solute = MW_solute_g_mol
        self.MW_solvent = MW_solvent_g_mol
        self.Vs = Vs
        self.vantHoffi=vantHoffi
        
        self.stage_df_list = []
        
        self.num_stages = [1,]
       
        self.df = pd.DataFrame(data=np.nan*np.ndarray((1,len(template))), columns = template, index = list(('S'+str(i) for i in range(1))))
        
        self.df = self.df.iloc[np.where(self.df.index.isnull()==False)]
        # mb = ['Feed Flow','Feed TDS','Permeate Flow','Permeate TDS','Concentrate Flow', 'Concentrate TDS', 'Selectivity','Recovery']
        # module_stats = ['Membrane Active Area','Modules Per Vessel','Membrane Area Per Module','# of Modules','Modules Per Vessel','NumVessels','Annual Cost of Modules']
        # costs = ['High Pressure Pumps', 'Circulation Pumps', 'Cost of Vessels', 'Peripheral Hardware','Total System Cost',]
        
        
        # pdb.set_trace()
        
        density = 1  ## kg/L
        
        self.input_tds = input_tds
        self.input_flow = feed_rate_mt_year*1000/density/3.8/(350*1440)  ### GPM
        
        
        self.df.loc['S0','Feed Flow']=self.input_flow
        self.df.loc['S0','Feed TDS']=self.input_tds
        
       
        
        self.df.loc[:,'Selectivity']=self.selectivity
        
        self.df['Flux']=flux.min()
        

        self.optimization_df= pd.DataFrame(columns=['residual'], index =[])
        if self.material=='SS':
            
            
            self.df_system_sizing = pd.Series(index = ['Installation Factor', 'Membrane Vessel','Pumps',
            'Piping',
            'Valves + Instruments',
            'Tanks',
            'Pre Filtration',
            'Temperature Control', 'Total'], data = [1.6,0,0,1.23,0.54,0.32,0.67,0.21,2.97],  )
        elif self.material =='Hastelloy':
            
            
            self.df_system_sizing = pd.Series(index = ['Installation Factor', 'Membrane Vessel','Pumps',
            'Piping',
            'Valves + Instruments',
            'Tanks',
            'Pre Filtration',
            'Temperature Control', 'Total'], data = [0.4,5,0,1.23*6,0.54*6,0.32*2,0.67,0.21,17.54],  )
        else:
            print("Material Not Recognized")
            return
        

 
        
        if autooptimize: 
            self.optimize_self()
        else:
            self.num_stages = num_stages
            
            self.massbalance()
            self.calculate_membrane_area()
            self.calculate_energy()
            self.calculate_system_cost() 
        return
    
    

        
    def calc_solute_recovery(self):
        c=self.get_concentrate_output()
        f=self.get_feed_input()
        return c[0]*c[1]/(f[0]*f[1])
        
    def calc_solvent_recovery(self):
        p=self.get_permeate_output()
        f=self.get_feed_input()
        return p[0]*(1-p[1])/(f[0]*(1-f[1]))
    
    def met_solvent_recovery(self):
        if self.calc_solvent_recovery()>self.solvent_recovery_target:
            return True
        else:
            return False
        
    def met_solute_recovery(self):
        if self.calc_solute_recovery()>self.solute_recovery_target:
            return True
        else:
            return False
        
    def met_concentrate_tds(self):
        if self.get_concentrate_output()[1]>self.concentrate_tds_target:
            return True
        else:
            return False
        
    def met_permeate_tds(self):
        if self.get_permeate_output()[1]<self.permeate_tds_target:
            return True
        else:
            return False
        
    
        
    def optimize_self(self, speedy=False, 
                      max_allowable_passes_before_termination = 20,
                      max_allowable_stages_before_termination = 20):
        
        self.num_stages=[5,]
        number_of_stages=self.num_stages[0]
        
        try:
            self.stage_df['Feed Flow']
            
        except:
            print('reseting feed flow')
            self.stage_df = pd.DataFrame(columns = ['Feed Flow','Feed TDS','Permeate Flow', 'Permeate TDS',
                                                'Concentrate Flow','Concentrate TDS',
                                                 'OsmoticPressure_Pa', ],
                                index = np.arange(number_of_stages), dtype=float)

        self.stage_df.loc[0,'Feed Flow']=1
        self.stage_df.loc[0,'Feed TDS']=0.000001
        
        self.massbalance()

        print("Beginning Pass 1 build.  Will add stages until reaching concentrate tds target")
        
        
        self.df.loc['S0', 'Feed Flow'] = self.get_feed_input()[0] 
        self.df.loc['S0', 'Feed TDS'] = self.permeate_tds_target
        
        
         
                
        print("Build complete", self.num_stages)



           

        
        self.break_down_stage_df_to_pass_df()
        self.calculate_membrane_area()
        self.calculate_energy()
        self.calculate_system_cost() 
        
        return 
        
        

        
        
    def add_stage_in_front(self):
        
        self.num_stages[0]+=1
        i=min(self.stage_df.index)
        n=pd.DataFrame(columns = self.stage_df.columns, index=[i])
        n.loc[i,'Feed Flow'] = self.stage_df['Permeate Flow'][0]
        n.loc[i,'Feed TDS'] = self.stage_df['Permeate TDS'][0]
        n.loc[i,'Permeate Connections']=-1
        self.stage_df.loc[:,'Permeate Connections']+=1
        self.stage_df.loc[0,'Permeate Connections']=0
        
        self.stage_df.index+=1
        self.stage_df=pd.concat((n, self.stage_df))
        self.feed_input_stage+=1
        stage_i=0
        
        bypass_ratio = self.stage_df.loc[stage_i,'Feed TDS']*0.1
        
        
        R=8.314
        T=298
        
        (feed_flow, feed_tds, true_permeate_flow,
        true_permeate_tds, permeate_flow,
        permeate_tds, concentrate_flow,
        concentrate_tds) = solve_stage_feed_mixed_with_permeate(
                                                            self.max_tmp_Pa*0.9,
                                                            self.stage_df.loc[stage_i,'Feed Flow'],
                                                           self.stage_df.loc[stage_i,'Feed TDS'],
                                                            self.selectivity,self.Vs, self.vantHoffi,T,
                                                            self.MW_solute, self.MW_solvent,
                                                            bypass_ratio, 
                                                            )
            

        self.stage_df.loc[stage_i, 'Permeate Ratio'] = true_permeate_flow/feed_flow

        self.stage_df.loc[stage_i,'Permeate TDS']=permeate_tds
        self.stage_df.loc[stage_i,'Permeate Flow'] = permeate_flow
        self.stage_df.loc[stage_i,'Concentrate Flow']=concentrate_flow
        self.stage_df.loc[stage_i,'Concentrate TDS'] = concentrate_tds
       
        self.stage_df.loc[stage_i,'True Permeate TDS']=true_permeate_tds
        self.stage_df.loc[stage_i,'True Permeate Flow'] =true_permeate_flow
        self.stage_df.loc[stage_i,'True Permeate Ratio'] =true_permeate_flow/self.stage_df.loc[stage_i,'Feed Flow']
        self.stage_df.loc[stage_i,'Bypass Ratio'] = (permeate_flow-true_permeate_flow)/(self.stage_df.loc[stage_i,'Feed Flow'])
        
    def add_stage_in_back(self):
        """adds a stage to the end of the stage_df.  Increases the final concentrate tds"""
        
        self.num_stages[0]+=1
        i=max(self.stage_df.index)
        
        self.stage_df.loc[i+1,'Feed Flow'] = self.stage_df['Permeate Flow'][i]
        self.stage_df.loc[i+1,'Feed TDS'] = self.stage_df['Permeate TDS'][i]
        self.stage_df.loc[i+1,'Permeate Connections']=i
        
        
        
        
        
        stage_i=i
        
        bypass_ratio = self.stage_df.loc[i,'Feed TDS']*0.1
        
        
        R=8.314
        T=298
        
        (feed_flow, feed_tds, true_permeate_flow,
        true_permeate_tds, permeate_flow,
        permeate_tds, concentrate_flow,
        concentrate_tds) = solve_stage_feed_mixed_with_permeate(
                                                            self.max_tmp_Pa*0.9,
                                                            self.stage_df.loc[stage_i,'Feed Flow'],
                                                           self.stage_df.loc[stage_i,'Feed TDS'],
                                                            self.selectivity,self.Vs, self.vantHoffi,T,
                                                            self.MW_solute, self.MW_solvent,
                                                            bypass_ratio, 
                                                            )
            

        self.stage_df.loc[stage_i, 'Permeate Ratio'] = true_permeate_flow/feed_flow

        self.stage_df.loc[stage_i,'Permeate TDS']=permeate_tds
        self.stage_df.loc[stage_i,'Permeate Flow'] = permeate_flow
        self.stage_df.loc[stage_i,'Concentrate Flow']=concentrate_flow
        self.stage_df.loc[stage_i,'Concentrate TDS'] = concentrate_tds
       
        self.stage_df.loc[stage_i,'True Permeate TDS']=true_permeate_tds
        self.stage_df.loc[stage_i,'True Permeate Flow'] =true_permeate_flow
        self.stage_df.loc[stage_i,'True Permeate Ratio'] =true_permeate_flow/self.stage_df.loc[stage_i,'Feed Flow']
        self.stage_df.loc[stage_i,'Bypass Ratio'] = (permeate_flow-true_permeate_flow)/(self.stage_df.loc[stage_i,'Feed Flow'])
        
        return
        
        
        
        
    
   
    def consolidate_permeate_connectivity(self):
        """consolidates the many passese into just a few passes
        num_passes gives the target number of passes"""
        
        num_passes = 6
        
        if len(self.stage_df['Permeate TDS'])<10:
            return
        
        cutoff_concentrations = [self.stage_df['Permeate TDS'][0]]
        cutoff_locations=[0]
         
        for target in np.linspace(min(self.stage_df['Feed TDS']), max(self.stage_df['Feed TDS']), num_passes)[0:num_passes-1]:
        # for target in np.linspace(min(self.stage_df['Permeate TDS']), max(self.stage_df['Permeate TDS']), 5)[0:4]:
            z=np.argmax(self.stage_df['Permeate TDS']>target)-1
            if z==-1:
                z=0
            cutoff_concentrations.append(self.stage_df['Feed TDS'][z])
            cutoff_locations.append(z)

        for c, l in zip(cutoff_concentrations, cutoff_locations):
            # pdb.set_trace()
            self.stage_df.loc[self.stage_df['Permeate TDS']>c, 'Permeate Connections']=l
        
        
        
    def set_permeate_connectivity(self):
        if 'Permeate Connections' in self.stage_df.columns:
            initial_permeate_connections = self.stage_df['Permeate Connections'].to_numpy()
        permeate_connections =np.array([0,])
        for i in range(1,len(self.stage_df.index)):
            p_tds = self.stage_df['Permeate TDS'].iloc[i]
            
            # permeate_cnct_stage =np.where((self.stage_df['Feed TDS'].iloc[1:]>p_tds).to_numpy() & (self.stage_df['Feed TDS'].iloc[:-1]<p_tds).to_numpy())[0]
            
            permeate_cnct_stage=np.intersect1d(np.where(self.stage_df['Feed TDS'].iloc[1:]>p_tds)[0],
                                        np.where(self.stage_df['Feed TDS'].iloc[:-1]<p_tds)[0])
            
            if len(permeate_cnct_stage)==0:
                permeate_cnct_stage=0
            else:
                permeate_cnct_stage=permeate_cnct_stage[0]
            permeate_connections=np.append(permeate_connections,permeate_cnct_stage)
           
        self.stage_df['Permeate Connections'] = permeate_connections
        # if np.all(initial_permeate_connections ==  permeate_connections):
        #     print("No changes in permeate connections")
        # else:
        #     print('Permeate Connections changed.')
        
        
    def massbalancestage(self, pass_num,number_of_stages):

        """compute mass balances for a permeate purification system
        
        In version July15 and beyond, this sets up a counter flow stage.  First the mass balances are solved assuming a tds of 0 a thte permeate,
        then the system iterates and solves for the values of recover that give self-consistent flows and concentrations """
        
        
        for currentstage in np.arange(number_of_stages):
            
            
            ##### set up the stage df
            stage_i = currentstage
            stage_iminus1 = currentstage-1
            try:
                # if currentstage==10:
                #     pdb.set_trace()
                if len(self.stage_df.index)<currentstage+1:
                    self.stage_df.loc[stage_i,'Feed TDS']=self.stage_df.loc[stage_iminus1,'Concentrate TDS']
                    self.stage_df.loc[stage_i,'Concentrate TDS']=self.stage_df.loc[stage_i,'Feed TDS']
                    self.stage_df.loc[stage_i,'Permeate TDS']=self.stage_df['Feed TDS'][stage_i]-self.stage_df['Concentrate TDS'][stage_i]
                    self.stage_df.loc[stage_i,'Concentrate Flow']=0.5
                    self.stage_df.loc[stage_i,'Permeate Flow']=0.5
                    self.stage_df.loc[stage_i,'Feed Flow']=1
                    
                if np.isnan(self.stage_df.loc[stage_i,'Concentrate TDS']):
                    self.stage_df.loc[stage_i,'Concentrate TDS']=self.stage_df['Feed TDS'][stage_i]
                    self.stage_df.loc[stage_i,'Concentrate Flow']=self.stage_df['Feed Flow'][stage_i]*0.5
                if np.isnan(self.stage_df['Permeate TDS'][stage_i]):
                    self.stage_df.loc[stage_i,'Permeate Flow']=self.stage_df['Feed Flow'][stage_i]*0.5
                    self.stage_df.loc[stage_i,'Permeate TDS']=self.stage_df['Feed TDS'][stage_i]
            except:
                raise
            
            
            ### set feed flow
            if currentstage ==0:
                
                self.stage_df.loc[0,'Feed Flow']=self.df['Feed Flow'][pass_num]
                self.stage_df.loc[0,'Feed TDS']=self.df['Feed TDS'][pass_num]
                
            elif currentstage>0:
    
                self.stage_df.loc[stage_i,'Feed Flow']=self.stage_df['Concentrate Flow'].iloc[stage_iminus1]
                self.stage_df.loc[stage_i,'Feed TDS']=self.stage_df['Concentrate TDS'].iloc[stage_iminus1]
                
                
            
            """solve current stage"""
           

            bypass_ratio = self.stage_df.loc[stage_i,'Feed TDS']*0.1
            
            
            R=8.314
            T=298
            
            (feed_flow, feed_tds, true_permeate_flow,
            true_permeate_tds, permeate_flow,
            permeate_tds, concentrate_flow,
            concentrate_tds) = solve_stage_feed_mixed_with_permeate(
                                                                self.max_tmp_Pa*0.9,
                                                                self.stage_df.loc[stage_i,'Feed Flow'],
                                                               self.stage_df.loc[stage_i,'Feed TDS'],
                                                                self.selectivity,self.Vs, self.vantHoffi,T,
                                                                self.MW_solute, self.MW_solvent,
                                                                bypass_ratio, 
                                                                )
                
  
            self.stage_df.loc[stage_i, 'Permeate Ratio'] = true_permeate_flow/feed_flow

            self.stage_df.loc[stage_i,'Permeate TDS']=permeate_tds
            self.stage_df.loc[stage_i,'Permeate Flow'] = permeate_flow
            self.stage_df.loc[stage_i,'Concentrate Flow']=concentrate_flow
            self.stage_df.loc[stage_i,'Concentrate TDS'] = concentrate_tds
           
            self.stage_df.loc[stage_i,'True Permeate TDS']=true_permeate_tds
            self.stage_df.loc[stage_i,'True Permeate Flow'] =true_permeate_flow
            self.stage_df.loc[stage_i,'True Permeate Ratio'] =true_permeate_flow/self.stage_df.loc[stage_i,'Feed Flow']
            self.stage_df.loc[stage_i,'Bypass Ratio'] = (permeate_flow-true_permeate_flow)/(self.stage_df.loc[stage_i,'Feed Flow'])
        
        
     
        
        
        for stage_i in np.arange(number_of_stages, ):   
            osmotic_pressure = calc_osmotic_pressure(self.stage_df['Concentrate TDS'][stage_i], self.stage_df['Permeate TDS'][stage_i],
                                      self.MW_solute, self.MW_solvent, i = self.vantHoffi,Vs=self.Vs, T=T )
            self.stage_df.loc[stage_i,'OsmoticPressure_Pa']=osmotic_pressure.round(0)
   
    def massbalance(self,setfeedflow=False, setfeedtds=False, passback_config=True):
        """compute mass balances for a permeate purification system"""
        
        print('balancing mass')
        residual_feedflow=10
        alpha=1
        if setfeedflow:
            self.df.loc['S0','Feed Flow']=self.get_feed_input()[0]
        
        
        if len(self.df.index)<len(self.num_stages):
            pass_name = 'S'+str(len(self.num_stages)-1)
            previous_pass = self.df.index[-1]
            self.df.loc[pass_name,'Feed Flow']=self.df['Permeate Flow'][previous_pass]
            self.df.loc[pass_name,'Feed TDS']=self.df['Permeate TDS'][previous_pass]
  
            
        ###At the top of this for-loop we erase the stage_df_list and rebuild it.   

        
        for currentpass in range(len(self.num_stages)):
            pass_name = 'S'+str(currentpass)
            previous_pass = 'S'+str(currentpass-1) 
            next_pass = 'S'+str(currentpass+1)
            
            self.massbalancestage(currentpass,self.num_stages[currentpass])
            self.stage_df_list.append(self.stage_df.copy()) 
            print(self.num_stages)
        

        cycle=1
        self.feed_input_stage=2
        self.set_permeate_connectivity()
            
        permeate_flows=[]
        

        
        while True:
            cycle+=1
            print(cycle)
            self.set_permeate_connectivity()
            
            self.set_feed_flows_and_conc_from_perm_and_conc(alpha=0.2)
            
            
            
            for stage_i in range(self.num_stages[0]):
                
                bypass_ratio = self.stage_df.loc[stage_i,'Feed TDS']*0.1
                R=8.314
                T=298
                (feed_flow, feed_tds, true_permeate_flow,
                true_permeate_tds, permeate_flow,
                permeate_tds, concentrate_flow,
                concentrate_tds) = solve_stage_feed_mixed_with_permeate(
                                                                    self.max_tmp_Pa*0.9,
                                                                    self.stage_df.loc[stage_i,'Feed Flow'],
                                                                    self.stage_df.loc[stage_i,'Feed TDS'],
                                                                    self.selectivity,self.Vs, self.vantHoffi,T,
                                                                    self.MW_solute, self.MW_solvent,
                                                                    bypass_ratio, 
                                                                    )

                self.stage_df.loc[stage_i, 'Permeate Ratio'] = true_permeate_flow/feed_flow

                self.stage_df.loc[stage_i,'Permeate TDS']=permeate_tds
                self.stage_df.loc[stage_i,'Permeate Flow'] = permeate_flow
                self.stage_df.loc[stage_i,'Concentrate Flow']=concentrate_flow
                self.stage_df.loc[stage_i,'Concentrate TDS'] = concentrate_tds
               
                self.stage_df.loc[stage_i,'True Permeate TDS']=true_permeate_tds
                self.stage_df.loc[stage_i,'True Permeate Flow'] =true_permeate_flow
                self.stage_df.loc[stage_i,'True Permeate Ratio'] =true_permeate_flow/self.stage_df.loc[stage_i,'Feed Flow']
                self.stage_df.loc[stage_i,'Bypass Ratio'] = (permeate_flow-true_permeate_flow)/(self.stage_df.loc[stage_i,'Feed Flow'])
            
           
        
            ###################################################################################################################
            f,p,c, permeate_connections, feed_input_stage=solve_flows(self.stage_df, self.get_feed_input()[0],
                                                    self.get_feed_input()[1],
                                                    feed_input_stage=None)
            
            
            self.feed_input_stage=feed_input_stage
            
            
            convergence_metric1 = f-self.stage_df['Feed Flow'].to_numpy()
            
            if np.all(convergence_metric1<0.0005) and self.set_feed_flows_and_conc_from_perm_and_conc(check_mode=True):
                break
            
            
            alpha = 0.1 ## high number for high damping
            self.stage_df['Feed Flow']=f*(1-alpha)+alpha*self.stage_df['Feed Flow']
            self.stage_df['Permeate Flow']=p*(1-alpha)+alpha*self.stage_df['Permeate Flow']
            self.stage_df['Concentrate Flow']=self.stage_df['Feed Flow']-self.stage_df['Permeate Flow']
            
            permeate_flows.append(self.stage_df['Permeate Flow'].to_numpy())
            
        
        
        ###############################################################################
        
    
    
        # return
        
        
      
        while True:
            print('dd')
            if not self.met_permeate_tds():
                 self.add_stage_in_front()
            if self.stage_df['Permeate TDS'].to_numpy()[1]<self.permeate_tds_target:
                self.stage_df = self.stage_df.drop(index = 0)
                self.stage_df.index=np.arange(len(self.stage_df))
                print('dropping stage',0 )
                self.num_stages[0]-=1
                self.feed_input_stage-=1
                self.stage_df['Permeate Connections']-=1
            if self.stage_df['Concentrate TDS'].to_numpy()[-2]>self.concentrate_tds_target:
                self.stage_df = self.stage_df.drop(index = self.stage_df.index[-1])
                print('dropping stage',self.stage_df.index[-1] )
                self.num_stages[0]-=1
                self.feed_input_stage-=1
                self.stage_df['Permeate Connections']-=1
            if not self.met_concentrate_tds():
                self.add_stage_in_back()
            
            
            
        
            #### now rerun but dont reset perm connectivity
            convergence_metric1=1
            cycle=1
            f_old=1
            while True:
                cycle+=1
                if cycle%10==0:
                    #
                    print(cycle, end='\r')
                
                
                print(self.stage_df['Permeate Connections'].to_numpy())
                
                
                self.set_feed_flows_and_conc_from_perm_and_conc(alpha=0.2)
                convergence_metric1 =  np.abs(f_old-self.stage_df['Feed Flow'].to_numpy())
                
                if np.all(convergence_metric1<0.0001 ) and self.set_feed_flows_and_conc_from_perm_and_conc(check_mode=True) :
                    
                    break
                f_old = self.stage_df['Feed Flow'].to_numpy()*1
                
                for stage_i in range(self.num_stages[0]):
                    
                    bypass_ratio = self.stage_df.loc[stage_i,'Feed TDS']*0.1
                    R=8.314
                    T=298
                    (feed_flow, feed_tds, true_permeate_flow,
                    true_permeate_tds, permeate_flow,
                    permeate_tds, concentrate_flow,
                    concentrate_tds) = solve_stage_feed_mixed_with_permeate(
                                                                        self.max_tmp_Pa*0.9,
                                                                        self.stage_df.loc[stage_i,'Feed Flow'],
                                                                       self.stage_df.loc[stage_i,'Feed TDS'],
                                                                        self.selectivity,self.Vs, self.vantHoffi,T,
                                                                        self.MW_solute, self.MW_solvent,
                                                                        bypass_ratio, 
                                                                        )
    
                    self.stage_df.loc[stage_i, 'Permeate Ratio'] = true_permeate_flow/feed_flow
    
                    self.stage_df.loc[stage_i,'Permeate TDS']=permeate_tds
                    self.stage_df.loc[stage_i,'Permeate Flow'] = permeate_flow
                    self.stage_df.loc[stage_i,'Concentrate Flow']=concentrate_flow
                    self.stage_df.loc[stage_i,'Concentrate TDS'] = concentrate_tds
                   
                    self.stage_df.loc[stage_i,'True Permeate TDS']=true_permeate_tds
                    self.stage_df.loc[stage_i,'True Permeate Flow'] =true_permeate_flow
                    self.stage_df.loc[stage_i,'True Permeate Ratio'] =true_permeate_flow/self.stage_df.loc[stage_i,'Feed Flow']
                    self.stage_df.loc[stage_i,'Bypass Ratio'] = (permeate_flow-true_permeate_flow)/(self.stage_df.loc[stage_i,'Feed Flow'])
                    
                    osm=osmotic_pressure = calc_osmotic_pressure(self.stage_df['Concentrate TDS'][stage_i], self.stage_df['Permeate TDS'][stage_i],
                                              self.MW_solute, self.MW_solvent, i = self.vantHoffi,Vs=self.Vs, T=T )
                    
            
                ###################################################################################################################
                
                
                
                f,p,c, permeate_connections, feed_input_stage=solve_flows(self.stage_df, self.get_feed_input()[0],
                                                        self.get_feed_input()[1],
                                                        feed_input_stage=self.feed_input_stage)
                
                
                
                
                
                alpha = 0.2 ## high number for high damping
                self.stage_df['Feed Flow']=f*(1-alpha)+alpha*self.stage_df['Feed Flow']
                self.stage_df['Permeate Flow']=p*(1-alpha)+alpha*self.stage_df['Permeate Flow']
                self.stage_df['Concentrate Flow']=self.stage_df['Feed Flow']-self.stage_df['Permeate Flow']
                
                permeate_flows.append(self.stage_df['Permeate Flow'].to_numpy())   
                
                if cycle>1000:
                    print('exceeded cycle limits')
                    break
            
            
                ###################################################################################################################
                
     
            if self.met_concentrate_tds() and self.met_permeate_tds(): 
                break      
        for stage_i in np.arange(len(self.stage_df)):   
            osmotic_pressure = calc_osmotic_pressure(self.stage_df['Concentrate TDS'][stage_i], self.stage_df['Permeate TDS'][stage_i],
                                      self.MW_solute, self.MW_solvent, i = self.vantHoffi,Vs=self.Vs, T=T )
            self.stage_df.loc[stage_i,'OsmoticPressure_Pa']=osmotic_pressure.round(0)
         
            
        # self.df.loc['S0','Feed Flow']=self.stage_df['Feed Flow'].iloc[0]
        # self.df.loc['S0','Permeate Flow']=self.stage_df['Permeate Flow'].iloc[0]
        # self.df.loc['S0','Concentrate Flow']=self.stage_df['Concentrate Flow'].iloc[-1]
        
        # self.df.loc['S0','Feed TDS']=self.stage_df['Feed TDS'].iloc[0]
        # self.df.loc['S0','Permeate TDS']=self.stage_df['Permeate TDS'].iloc[0]
        # self.df.loc['S0','Concentrate TDS']=self.stage_df['Concentrate TDS'].iloc[-1]
            
            
        
    def save_configuration(self, filename):
        
        # Creating Excel Writer Object from Pandas  
        writer = pd.ExcelWriter(filename,engine='xlsxwriter')   
        workbook=writer.book
        
        for i,df in enumerate(self.stage_df_list):
            
          
            df.to_excel(writer,sheet_name='Pass'+str(i),startrow=0 , startcol=0)   
        writer.close()
            
        
    def get_system_structure(self):
        return '-'.join(list((str(i) for i in self.num_stages)))
    def get_number_of_passes(self):
        return self.df.shape[0]
    
    
    def set_feed_flows_and_conc_from_perm_and_conc(self, alpha=0.2,check_mode=False, verbose=False):
        """sums the flows and calculated the concentrations of the combined permeate flows for connected stages"""
        
        
        self.stage_df.loc[0,'Permeate Connections'] = -1  
        new_feed_flow =[]
        new_feed_tds =[]
        for i in self.stage_df.index:
            if i>0:
                prior_stage_conc_flow = self.stage_df['Concentrate Flow'][i-1] 
                prior_stage_conc_tds = self.stage_df['Concentrate TDS'][i-1]
            else:
                prior_stage_conc_flow=0
                prior_stage_conc_tds=0.5
                
            perms = self.stage_df.loc[self.stage_df['Permeate Connections']==i]
            if len(perms)==0:
                perm_flow=0
                perm_tds=0.5
            else:
                perm_flow = perms['Permeate Flow'].sum()
                if perm_flow==0:
                    perm_tds=0
                else:
                    perm_tds = (perms['Permeate Flow']*perms['Permeate TDS']).sum()/perm_flow
            
            if self.feed_input_stage==i:
                feed_input_flow = self.get_feed_input()[0]
                feed_input_tds = self.get_feed_input()[1]
            else:
                feed_input_flow=0
                feed_input_tds=0.5
            new_feed_flow.append( perm_flow+prior_stage_conc_flow+feed_input_flow)
            new_feed_tds.append(np.sum(perm_tds*perm_flow + prior_stage_conc_flow*prior_stage_conc_tds+feed_input_flow*feed_input_tds)/new_feed_flow[-1])
            
        
        if check_mode:
            
            if np.all(np.abs(self.stage_df['Feed Flow']-new_feed_flow)<0.005) and  np.all(np.abs(self.stage_df['Feed TDS']-new_feed_tds)<0.0005):
                return True
            else:
                if verbose:
                    print("new Feed Flow")
                    print(new_feed_flow)
                    print('new feed tds')
                    print(new_feed_tds)
                return False
            
       
        
        self.stage_df.loc[:,'Feed Flow']=  alpha*np.array(new_feed_flow) + (1-alpha)*self.stage_df.loc[:,'Feed Flow']
        self.stage_df.loc[:,'Feed TDS']=   alpha*np.array(new_feed_tds) + (1-alpha)*self.stage_df.loc[:,'Feed TDS']
        return
        
        return 
    
    def break_down_stage_df_to_pass_df(self):
        """takes the single long stage_df with permeate connections and creates both stage_df_list and df"""
        self.stage_df['Master Index']=self.stage_df.index
        self.stage_df_list=[]        
        self.stage_df.loc[0,'Permeate Connections'] = -1  
        unique_passes = self.stage_df['Permeate Connections'].unique()
        unique_passes=np.append(unique_passes,max(self.stage_df.index)+1)
        pass_num=0
        for i in range(len(unique_passes)-1):
            ### note, don't group the passes by their permeate connectoins.  Instead, a pass begins at each unique permeate connection.
            self.stage_df.loc[unique_passes[i]:unique_passes[i+1]-1, 'Pass']=pass_num
            pass_num+=1
            if i>0:
                self.stage_df_list.append(self.stage_df.loc[unique_passes[i]:unique_passes[i+1]-1].reset_index())
        
        self.stage_df_list.reverse()
        self.df = self.df.reset_index()
        self.df=self.df.drop(labels='index', axis=1)
        
        
        for i in range(len(self.stage_df_list)):
            print('df ',i)
            self.df.loc[i,['Feed Flow', 'Feed TDS']]=self.stage_df_list[i][['Feed Flow','Feed TDS']].iloc[0]
            self.df.loc[i,['Concentrate Flow', 'Concentrate TDS']]=self.stage_df_list[i][['Concentrate Flow','Concentrate TDS']].iloc[-1]
            self.df.loc[i,'Permeate TDS']=(self.stage_df_list[i]['Permeate Flow']*self.stage_df_list[i]['Permeate TDS']).sum()/self.stage_df_list[i]['Permeate Flow'].sum()
            self.df.loc[i,'Permeate Flow']=self.stage_df_list[i]['Permeate Flow'].sum()
            self.df.loc[i,'MaxOsmPressure_psi'] = self.stage_df_list[i]['OsmoticPressure_Pa'].max()/101325*14.7
            
        self.num_stages = list((len(s) for s in self.stage_df_list))
        # pdb.set_trace()
        
    def calculate_membrane_area(self):
        
            
            
        for s in self.stage_df_list:
            s['Flux'] = self.df.loc[0, 'Flux']
            s['Modules Per Vessel']=6
            s['Membrane Area Per Module']=230
            
            s['Module Lifetime']=1 # year
            s['Unit Cost of Modules'] = 1500
            s['Unit Cost of Vessels'] = 5343
            s['Circ Pump Efficiency']=0.7
            s['HP Pump Efficiency']=0.7
            s['Operating Pressure']= self.max_tmp_Pa/101325*14.7  ## psi
            s['Pressure Drop']=5 ## psi
            s['Circ Flow Per Vessel'] = 5 ## GPM

            z=s.columns
            
            s['Membrane Active Area'] = s['True Permeate Flow'].to_numpy()/(s['Flux']/24/60)
            s['# of Modules']  = (s['Membrane Active Area'].to_numpy()/s['Membrane Area Per Module']).apply(np.ceil)
            s['NumVessels'] = (s['# of Modules'].to_numpy()/s['Modules Per Vessel']).apply(np.ceil)
            s['Annual Cost of Modules'] = s['# of Modules'].to_numpy()*s['Unit Cost of Modules']/s['Module Lifetime']
        
    def calculate_system_cost(self):

            
        for s in self.stage_df_list:
            
            
            mini_vessel = np.where(s['# of Modules'].to_numpy()==1)[0]
            s.loc[mini_vessel, 'Unit Cost of Vessels']/=4
            s['Cost of Vessels']= s['Unit Cost of Vessels'].to_numpy() *  s['NumVessels']
            
            
            
            s.loc[:, 'High Pressure Pumps']=0
            add_hp_pump =np.where(s['HP Pump Power'].to_numpy()>0.01)[0]
            # s.loc[add_hp_pump, 'High Pressure Pumps']=(s.loc[add_hp_pump,'HP Pump Power']*391.78+54982) * 1.12  ### factor of 12% added for auxilliary pumps
            s.loc[add_hp_pump, 'High Pressure Pumps']=(s.loc[add_hp_pump,'HP Pump Power']*1400+8000) * 1.12  ### factor of 12% added for auxilliary pumps
            
            s['Circulation Pumps']=0
            add_circ_pump =np.where(s['Circ Pump Power'].to_numpy()>0.01)[0]
            # s.loc[add_circ_pump, 'Circulation Pumps']=(s.loc[add_circ_pump,'Circ Pump Power']*170.97+47072)* 1.12  ### factor of 12% added for auxilliary pumps
            s.loc[add_circ_pump, 'Circulation Pumps']=(s.loc[add_circ_pump,'Circ Pump Power']*1400+8000)* 1.12  ### factor of 12% added for auxilliary pumps
            
            s['Counterflow Circulation Pumps']=0
            add_circ_pump =np.where(s['Counterflow Circ Pump Power'].to_numpy()>0.01)[0]
            s.loc[add_circ_pump, 'Counterflow Circulation Pumps'] = (s.loc[add_circ_pump, 'Counterflow Circ Pump Power']*1400+8000)* 1.12  ### factor of 12% added for auxilliary pumps
            
            s['Peripheral Hardware'] = self.df_system_sizing.loc['Total']*s['Cost of Vessels'].values
            
            
            s['Total System Cost'] = s['Cost of Vessels'].to_numpy()+s['High Pressure Pumps']+s['Circulation Pumps']+s['Peripheral Hardware']+s['Counterflow Circulation Pumps']
            s['Construction and Assembly']= self.df_system_sizing.loc['Installation Factor']*s['Total System Cost'].to_numpy()
            s['Total System Cost']*=(1+self.df_system_sizing['Installation Factor'])
        
        
    
        """feed rate"""
        GALtoM3=3.8/1000
        PSI_to_PA = 101325/14.7
        for ind, s in enumerate(self.stage_df_list):
            
            
            
            self.df.loc[ind,'HP Pump Power']=s['HP Pump Power'].sum()
            self.df.loc[ind,'Circ Flow Rate'] = s['Circ Flow Rate'].sum()
            self.df.loc[ind, 'Circ Pump Power']=s['Circ Pump Power'].sum()
            
            for col in ['NumVessels', 'Cost of Vessels',
                        'High Pressure Pumps',
                        'Circulation Pumps',
                        'Peripheral Hardware',
                        'Total System Cost',
                        'High Pressure Pumps',
                        'Construction and Assembly','Total System Cost','# of Modules','Circulation Pumps']:
                
                if col not in self.df.columns:
                    self.df[col]=np.nan
                self.df.loc[ind,col] =s[col].sum()
            
        # pdb.set_trace()   
        return
                                        
                                        
    def get_system_cost(self):
        
        return np.sum(self.df['Total System Cost'])
    
    def get_annual_module_cost(self):
        return sum(stage_df['Annual Cost of Modules'].sum() for stage_df in self.stage_df_list)
        # return np.sum(self.df['Annual Cost of Modules'])
    
    def get_levelized_cost(self, verbose = False):
        Electricity_Cost=0.018
        system_cost = self.get_system_cost()
        labor =80000
        annualmembranecost =  self.get_annual_module_cost()
        annualeleccosts = self.get_total_power()*350*24*Electricity_Cost
        if annualeleccosts==-np.inf:
            # print('annual elec costs is infinite')
            annualeleccosts=0
        if self.levelized_cost_basis =='concentrate':
            ProductRate = self.get_concentrate_output()[0]*350*24*60*3.8/1000  ## MT/year
        
        if self.levelized_cost_basis =='permeate':
            ProductRate = self.get_permeate_output()[0]*350*24*60*3.8/1000  ## MT/year
        
        if self.levelized_cost_basis =='feed':
            ProductRate = self.get_feed_input()[0]*350*24*60*3.8/1000  ## MT/year
       
        opex = annualmembranecost + labor + annualeleccosts
                         
        levelized_cost = (system_cost+opex*10)/ (ProductRate*10)
        if verbose:
            print("Levelized cost is ${0:.2f}/MT of {1}".format(levelized_cost, self.levelized_cost_basis))
        return levelized_cost
    
    def get_stats(self):
        mb = ['Feed Flow','Feed TDS','Permeate Flow','Permeate TDS','Concentrate Flow', 'Concentrate TDS', 'Selectivity',]
        module_stats = ['Membrane Active Area','Modules Per Vessel','Membrane Area Per Module','# of Modules','Modules Per Vessel','NumVessels','Annual Cost of Modules']
        costs = ['High Pressure Pumps', 'Circulation Pumps', 'Cost of Vessels', 'Peripheral Hardware','Total System Cost',]
        energy = ['HP Pump Power', 'Circ Pump Power',]
        print(self.df[mb])  
        print(self.df[module_stats])
        print(self.df[costs])
        print(self.df[energy])
        
    
    
    def get_total_power(self, unit='kW'):
        """returns power in kilowatts"""
        if unit =='kW':
            return self.df[['HP Pump Power', 'Circ Pump Power' ]].sum().sum()
        elif unit =='W':
            return self.df[['HP Pump Power', 'Circ Pump Power' ]].sum().sum()*1000
    
    def get_feed_input(self):
        return self.input_flow,self.input_tds
    
    
    def get_concentrate_output(self):
        
        return self.stage_df['Concentrate Flow'].iloc[-1],self.stage_df['Concentrate TDS'].iloc[-1]
        # return self.df['Concentrate Flow'].sum(), (self.df['Concentrate Flow']*self.df['Concentrate TDS']).sum()/self.df['Concentrate Flow'].sum()

        
    def get_permeate_output(self):
        
        return self.stage_df.iloc[0]['Permeate Flow'], self.stage_df.iloc[0]['Permeate TDS']
        #return self.df['Permeate Flow'][self.final_stage], self.df['Permeate TDS'][self.final_stage]
            
    def calculate_energy(self):
        """feed rate"""
        GALtoM3=3.8/1000
        PSI_to_PA = 101325/14.7
        for s in self.stage_df_list:
            

            s.loc[0,'HP Pump Power']=(s.loc[0,'Feed Flow']*GALtoM3*((s.loc[0,'Operating Pressure']+s.loc[0,'Pressure Drop'])*PSI_to_PA-101325)) \
                                                /(s.loc[0,'HP Pump Efficiency']*1000)/60
           
            s.loc[1:, 'HP Pump Power']=(s.loc[1:, 'Feed Flow']*GALtoM3*((s.loc[1:,'Pressure Drop'])*PSI_to_PA-101325)) \
                                                    /(s.loc[1:,'HP Pump Efficiency']*1000)/60
            
            
            s['Circ Flow Rate'] = s['NumVessels'].to_numpy()*s['Circ Flow Per Vessel']-s['Feed Flow']
            s.loc[s['Circ Flow Rate']<0,'Circ Flow Rate']=0
            s['Circ Pump Power']=((s['Circ Flow Rate'].to_numpy()*GALtoM3)*(s['Pressure Drop'].to_numpy()*PSI_to_PA))  \
                                            /(s['Circ Pump Efficiency'].to_numpy()*1000)/60
            
            s['Counterflow Circ Flow Rate']=  s['NumVessels'].to_numpy()*s['Circ Flow Per Vessel']-s['Permeate Flow']
            # print("Counterflow circ flow rate is artificially set to zero")
            s.loc[s['Counterflow Circ Flow Rate']<0,'Counterflow Circ Flow Rate']=0
            s['Counterflow Circ Pump Power']=((s['Counterflow Circ Flow Rate'].to_numpy()*GALtoM3)*(s['Pressure Drop'].to_numpy()*PSI_to_PA))  \
                                            /(s['Circ Pump Efficiency'].to_numpy()*1000)/60
                                            
            
                                        
    def check_solution(self):
        
        passing = True
        
        
        print("Checking for system overall mass balance")
        x= np.prod(self.get_feed_input())-np.prod(self.get_concentrate_output())-np.prod(self.get_permeate_output())
        p=np.prod(self.get_feed_input())
        if np.abs(x)> 0.01:
            print('.................Total System solute balance fail.  Error is {0:.4f}/{1:.4f}'.format( x,p))
        else:
            print('.................Total System solute balance PASS')
        x= self.get_feed_input()[0]-self.get_concentrate_output()[0]-self.get_permeate_output()[0]
        if np.abs(x)> 0.01:
            print('.................Total System flow balance fail. Error is  {0:.4f}/{1:.4f}'.format( x,self.get_feed_input()[0]))
        else:
            print('.................Total System flow balance PASS')
        
        print("Checking for correct mass balances in  stage_df")
        s =self.stage_df
            
        mass_balance_all = s['Permeate Flow'] + s['Concentrate Flow'] - s['Feed Flow']
        mass_balance_solute = s['Permeate Flow']*s['Permeate TDS'] + s['Concentrate Flow']*s['Concentrate TDS'] - s['Feed Flow']*s['Feed TDS']
        
        if np.any(mass_balance_all.abs()>0.1):
            print(".............mass_balance_all fail on stages", np.where(mass_balance_all.abs()>0.1)[0])
            passing=False
        else:
            print(".............mass_balance_all PASS")
        if np.any(mass_balance_solute.abs()>0.1):
            print(".............mass_balance_solute fail on stages",np.where(mass_balance_all.abs()>0.1)[0])
            passing=False
        else:
            print(".............mass_balance_solute PASS")
                
        print('Checking for mass balances in pass dataframe')        
        mass_balance_all = self.df['Permeate Flow'] + self.df['Concentrate Flow'] - self.df['Feed Flow']
        mass_balance_solute = self.df['Permeate Flow']*self.df['Permeate TDS'] + self.df['Concentrate Flow']*self.df['Concentrate TDS'] - self.df['Feed Flow']*self.df['Feed TDS']
        
        if np.any(mass_balance_all.abs()>0.1):
            print(".............mass_balance_all fail on self.df", np.where(mass_balance_all.abs()>0.1))
            passing=False
        else:
            print("............Did not detect any flow summation issues")
        if np.any(mass_balance_solute.abs()>0.1):
            print(".............mass_balance_solute fail on self.df", np.where(mass_balance_all.abs()>0.1))
            passing=False
        else:
            print("............Did not detect any solute summation issues")
        
        print("Confirming osmotic pressure does not exceed limits")
        try:
            osmotic_pressure = self.stage_df['OsmoticPressure_Pa']
            
            if np.any(osmotic_pressure> self.max_tmp_Pa):
                print('TMP fail in stage_df', np.where(osmotic_pressure> self.max_tmp_Pa))
                passing=False
            else:
                print("............No osmotic pressure issued dected")
        except:
            print("Could not calculate osmotic pressure")
            
        print("Confirming selectivity does not exceed limits")
        apparent_selectivity = s['Concentrate TDS']/s['Permeate TDS']
        error_arg = np.where(apparent_selectivity>self.selectivity)[0]
        if len(error_arg)>0:
            print('Selectivity violated in stage_df', error_arg)
            passing=False
        else:
            print("............PASS")
        
            
        print("Checking for agreement betweeen system data frame and pass dataframes ")
        local_pass=True
        for i in range(len(self.df.index)):
            
            if self.df.iloc[i]['Feed Flow'] != self.stage_df_list[i].iloc[0]['Feed Flow']:
                print("feed flow in df does not match feed flow in stage_df", i)
                passing=False
                local_pass=False
            if self.df.iloc[i]['Permeate Flow'] != self.stage_df_list[i]['Permeate Flow'].sum():
                print("Permeate flow in df does not match Permeate flow in stage_df", i)
                passing=False
                local_pass=False
            if self.df.iloc[i]['Concentrate Flow'] != self.stage_df_list[i].iloc[-1]['Concentrate Flow']:
                print("Concentrate flow in df does not match Concentrate flow in stage_df", i)
                passing=False
                local_pass=False
        if local_pass:
            print("....Agreement satisfactory.")
  
        
        
        if passing:
            print("...........Check complete.  Everything looks good.")
        else:
            print("...........errors detected.  Check for issues")
        return passing
                                        
                                        

    
def update_save_file(temp):
    
    df = pd.read_excel(save_file, index_col=0)
    df=pd.concat((df,temp))
    df.to_excel(save_file)
    
    
def create_data(feed_rate_mt_year):
    """feed_rate in MT/year"""
    
 
    
    input_tds = [0.019]#np.arange(0.2,1,0.1)
    concentrate_tds_target = [0.114]#0.7,0.8,0.9]
    solute_recovery_target =[0.60]
    solvent_recovery_target = [0]
    permeate_tds_target = [0.01]

    
    selectivity_range=[2]
    fluxes_GFD=[30]
    materials = ['SS']
    
    ### parameters for system tuning
    
    
    
    
    #### this block uses a single rejection and recovery value for each membrane
    rr = list(itertools.product(input_tds, 
                                solute_recovery_target,
                                solvent_recovery_target,
                                concentrate_tds_target,
                                permeate_tds_target,
                                selectivity_range,

                                fluxes_GFD,
                                materials,
                                
                                ))
    
    costs = pd.DataFrame( columns = ['input_tds',
                                     'solute_recovery_target',
                                     'solvent_recovery_target',
                                     'concentrate_tds_target',
                                     'permeate_tds_target',
                                     
                                     'sel1',
                                     'Flux','Material'], data = rr)
    costs['NumStagesPerPass']=""
   
    

    
    costs = costs.drop_duplicates()
    
    
    costs['SystemCost']=np.nan
    costs['PermeateTDS']=np.nan
    costs['Power']=np.nan
    print("WARNING!!! NEED TO UPDATE FOR ACCURATE PUMP COSTS!!!!!")
    Electricity_Cost = 0.18 ### $/kWhr
    costs.index = np.arange(len(costs.index))
    
    for i in costs.index:
        if i%100==0:
            print('index', i,'/', len(costs.index))
            
         
        try:
            v=Via(feed_rate_mt_year,
                costs.loc[i,'input_tds'],
                costs.loc[i,'solute_recovery_target'],
                costs.loc[i,'solvent_recovery_target'],
                costs.loc[i,'concentrate_tds_target'],
                costs.loc[i,'permeate_tds_target'],
                costs.loc[i,'sel1'],
                #costs.loc[i,['rec1','rec1','rec1','rec1']],
                  template = df_columns,
                  flux=costs.loc[i,'Flux'],
                  material=costs.loc[i,'Material'],
                  levelized_cost_basis='concentrate',
                  MW_solute_g_mol = 108,
                  MW_solvent_g_mol = 18,
                  rho_solvent_kg_m3 =1000,
                  rho_solute_kg_m3 =1000,
                  vantHoffi=2
                  )
        except:
            
            raise
        try:
            v.df=v.df.drop(labels=[np.nan,'Mass','Volume','Rejection','Tanks'], axis=1)
        except:
            pass
        
        return costs, v
        v.calculate_membrane_area()
        v.calculate_energy()
        v.calculate_system_cost()
        
        costs.loc[i,'SoluteRecovery'] = v.calc_solute_recovery()
        costs.loc[i,'SolventRecovery'] = v.calc_solvent_recovery()
        costs.loc[i,'NumVessels'] = v.df['NumVessels'].sum()
        costs.loc[i, 'NumPasses'] = v.get_number_of_passes()
        costs.loc[i,'NumStagesPerPass'] = '-'.join(str(i) for i in v.num_stages)
        
        costs.loc[i,'SystemCost'] = v.get_system_cost()
        costs.loc[i,'AnnualMembraneCost'] = v.get_annual_module_cost()
        costs.loc[i,'AnnualElecCosts']=v.get_total_power()*350*24*Electricity_Cost
        costs.loc[i, 'PermeateRate_MT_year'] = v.get_permeate_output()[0]*1915.2  ## MT/year
        costs.loc[i, 'PermeateTDS'] =v.get_permeate_output()[1]
        conc_output= v.get_concentrate_output()
        costs.loc[i, 'ConcentrateRate_MT_year']=conc_output[0]*1915.2  ## MT/year
        costs.loc[i, 'ConcentrateTDS']=conc_output[1]
        costs.loc[i,'ConcentrateTDS'] = costs.loc[i,'ConcentrateTDS'].round(5)
        costs.loc[i,'LCOSA']=v.get_levelized_cost()
        
        for col in ['Cost of Vessels','High Pressure Pumps','Circulation Pumps','Peripheral Hardware','Construction and Assembly']:
            costs.loc[i, col]=v.df[col].sum()
        
            

    # breakpoint
    costs = costs.loc[costs['NumVessels'].isnull()==False]
    
    costs['ConcentrateProcessingCosts_USD_year'] = costs['ConcentrateRate_MT_year']*200  ### dollars

    costs['LaborCosts'] = 18300
    costs['OPEXYearly_NoLabor'] =  costs['AnnualMembraneCost'] + costs['LaborCosts'] + costs['AnnualElecCosts']
    costs['CAPEXYearly'] = costs['SystemCost']/10
    costs['OPEXYearly'] =  costs['AnnualMembraneCost']+ costs['AnnualElecCosts']\
                                                + costs['LaborCosts']
                                                

    
    costs.to_excel("temp.xlsx")
    
    
    return costs, v

#

if __name__== "__main__":
    close('all')
    os.chdir('H:/Shared drives/Via Separations Main Folder/New Applications/SQI')
    feed_rate_mt_year =100000 ###MT/year of spent acid
    
    
    # template_df=pd.read_excel(template_file, sheet_name='VIA SYS 40 Remake',
    #               usecols = (2,3,4),skiprows = 5, nrows = 63 )
    # column_names = template_df['Unnamed: 2'].to_list()
    
    costs, v=create_data(feed_rate_mt_year)
    
    
    # v.save_configuration('H:/Shared drives/Via Separations Main Folder/New Applications/SQI/GreenMethanolSystemConfig.xlsx')
    import copy
    r=copy.deepcopy(v)
    r.df = copy.deepcopy(v.df)
    # v.num_stages=[15,]
    # v.df=v.df.iloc[0:len(v.num_stages)]
    # v.massbalance() 
    
    # v.calculate_membrane_area()
    # v.calculate_energy()
    # v.calculate_system_cost()
    
            
   
    
    
    """When permeate is the product, you will need total recovery and purity requirements.  """
     
