# -*- coding: utf-8 -*-
"""
Created on Thu Jun 27 13:23:05 2019

@author: rheil
"""
#==============================================================================
# Imports and assign paths
#==============================================================================
import sys
import os
import pandas as pd
import numpy as np

sys.path.append('D:/dev/glue-sb/')
import dirfuncs
dropbox_dir = dirfuncs.guess_dropbox_dir()
data_dir = dropbox_dir + 'soyM/analysis/3-11-20/'

pt_df = pd.read_csv(data_dir + 'wide.csv')

# =============================================================================
# Create folders for tables and figures
# =============================================================================
for dir_name in [data_dir + 'tables/', data_dir + 'figures/']:
    if not os.path.exists(dir_name):
        os.mkdir(dir_name)

# =============================================================================
# Convert to long table
# =============================================================================
pt_df = pt_df.loc[(pt_df['mb2_vfor_2000']==1)]
stubs = ['temp_', 'trmm_', 'gts_']
long_df = pd.wide_to_long(pt_df, stubnames = stubs, i = 'ptid', 
                          j = 'year')
                          
change_stubs = [stub for stub in stubs if '_'==stub[-1:]]
long_df = long_df.rename(columns = {col: col[:-1] for col in change_stubs})
long_df = long_df.reset_index()

long_df['mb2_for'] = 1
long_df.loc[long_df['mb2_for_2000']==0, 'mb2_for'] = np.nan
long_df.loc[long_df['year'] == long_df['mb2_y_defor'], 'mb2_for'] = 0
long_df.loc[long_df['year'] > long_df['mb2_y_defor'], 'mb2_for'] = np.nan
long_df['mb2_defor'] = long_df['mb2_for'].apply(lambda x: not(x)).astype(int)
long_df.loc[long_df['mb2_for'].isnull(), 'mb2_defor'] = np.nan

long_df['mb2_vfor'] = 1
long_df.loc[long_df['mb2_vfor_2000']==0, 'mb2_vfor'] = np.nan
long_df.loc[long_df['year'] == long_df['mb2_y_defor'], 'mb2_vfor'] = 0
long_df.loc[long_df['year'] > long_df['mb2_y_defor'], 'mb2_vfor'] = np.nan
long_df['mb2_vdefor'] = long_df['mb2_vfor'].apply(lambda x: not(x)).astype(int)
long_df.loc[long_df['mb2_vfor'].isnull(), 'mb2_vdefor'] = np.nan

# =============================================================================
# Assign policy treatments
# =============================================================================
#long_df.loc[long_df['biome']==1, 'legal_amazon'] = 1
long_df['soy_suit'] = (long_df['suit']>0).astype(int) * (long_df['GAEZsuit']>40).astype(int)

long_df['car_now'] = (long_df['year']>=long_df['car_year']).astype(int)
long_df['gts'] = long_df['gts'].fillna(0)
long_df['gts_now'] = ((long_df['year']>2007) * long_df['gts']).astype(int)
long_df['asm_now'] = ((long_df['year']>2005) * long_df['biome']).astype(int)

long_df.loc[long_df['year']==2007, 'gts'] = 0
long_df['gts_now'] = long_df['prodes_mon'] * long_df['gts'] * long_df['biome'] * \
    np.logical_not(long_df['set']).astype(int) * np.logical_not(long_df['pa']).astype(int)

#==============================================================================
# Save long data
#==============================================================================
long_df = long_df.sort_values('random')
long_df.to_csv(data_dir + 'long.csv', index = False)

# =============================================================================
# Create long dataset for soy conversion analysis (2 time period)
# =============================================================================
pt_df = pd.read_csv(data_dir + 'wide.csv')

## Soy suitability class
pt_df['soy_suit'] = (pt_df['suit']>0).astype(int) * (pt_df['GAEZsuit']>40).astype(int)

## Create starting soy identifiers
pt_df['a_start_soy_2006'] = pt_df['a_soy_2000']
pt_df['a_start_soy_2017'] = pt_df['a_soy_2006']

## Create starting forest identifiers
pt_df['mb_start_for_2006'] = pt_df['mb2_vfor_2000']
pt_df['nodefor_pre2006'] = ((pt_df['mb2_y_defor']>2006) | (pt_df['mb2_y_defor'].isnull())).astype(int)
pt_df['mb_start_for_2017'] = pt_df['mb2_vfor_2000'] * pt_df['nodefor_pre2006']

# Get mean temp and precip values over each period
for var in ['temp', 'trmm']:
    pt_df['m' + var + '_2006'] = pt_df[[var + '_' + str(y) for y in range(2000, 2007)]].astype(float).mean(axis = 1)   
    pt_df['m' + var + '_2017'] = pt_df[[var + '_' + str(y) for y in range(2007, 2017)]].astype(float).mean(axis = 1)   

# Convert to long dataset
stubs = ['a_soy_', 'mb_start_for_', 'a_start_soy_', 'mtemp_', 'mtrmm_']
long_df = pd.wide_to_long(pt_df, stubnames = stubs, i = 'ptid', 
                          j = 'year')
long_df = long_df.reset_index()
change_stubs = [stub for stub in stubs if '_'==stub[-1:]]
long_df = long_df.rename(columns = {col: col[:-1] for col in change_stubs})

# Export data for analysis in stata
long_df = long_df.sort_values('random')
out_df = long_df[['ptid', 'year', 'mb2_vfor_2000', 'mb2_y_defor', 'prodes_mon',
                  'propid', 'temp_2000', 'trmm_2000', 'urbandist', 'roaddist', 'soy_suit',
                  'municcode', 'state', 'a_soy', 'biome', 
                  'legal_amazon', 'set', 'random', 'dist_amb', 'dist_aml',
                  'pa', 'pa_noncar_elig', 'pa_car_elig', 'pa_indig', 'pa_pi', 'pa_us', 'pa_quil', 'pa_mil',
                  'mb_start_for', 'a_start_soy', 'mtemp', 'mtrmm', 'car_year']]
out_df.to_csv(data_dir + 'soy_conversion.csv', index = False)





# =============================================================================
# multinomial logit for on-property leakage
# =============================================================================
