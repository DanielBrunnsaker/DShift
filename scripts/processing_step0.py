#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Dec 11 11:00:25 2021

@author: danbru
"""
import pandas as pd
import os

def SNfilter(data,column_name,threshold):
    cleaned_data = data[pd.to_numeric(data[column_name]) > threshold]
    return(cleaned_data)

def Fillfilter(data,column_name,threshold):
    cleaned_data = data[pd.to_numeric(data[column_name]) > threshold]
    return(cleaned_data)

def blankfilter(data,regex_blanks, regex_samples):
    blanks_average = data.iloc[:,data.columns.str.contains(regex_blanks)].mean(axis = 1)
    sample_average = data.iloc[:,data.columns.str.contains(regex_samples)].mean(axis = 1)
    cleaned = data[sample_average.ge(blanks_average)]
    return(cleaned)

# ------ set parameters --------

sn_threshold = 10
fill_threshold = 0.1
threshold_0 = 0.66

# ------ set path --------

rel_path = '/Users/danbru/Desktop/DShift' #set project folder
os.chdir(rel_path)

# --------- load files, batch2 ---------
src_cols = pd.read_excel('data/metabolomics/Metabolomics_B2.xlsx').iloc[1:5,]
src_data = pd.read_excel('data/metabolomics/Metabolomics_B2.xlsx').iloc[4:,]
src_data.columns = src_data.iloc[0,]
src_data = src_data.iloc[1:,0:(src_data.shape[1]-2)]

# -------- filtering --------
column_name = 'S/N average'
data = SNfilter(src_data,column_name,sn_threshold)

column_name = 'Fill %'
data = Fillfilter(src_data,column_name,fill_threshold)

#filter out if mean is higher than blank
regex_samples = 'WT|DLD3|FAA1|GAL11|MEK1|OCA1|PCL1|RME|RTS3|TDA1|YGR'
regex_blanks = 'blank|MeOH'
data = blankfilter(data, regex_blanks, regex_samples)

#filter out zeros in each group
regex_b7 = 'B7|B5'
temp_dat = data.iloc[:,data.columns.str.contains(regex_b7)]
temp_dat = temp_dat.iloc[:,~(temp_dat.columns.str.contains('QC'))]
sumzeros = (temp_dat == 0).astype(int).sum(axis=1)
zeroBool_glucose = sumzeros.ge(threshold_0*len(temp_dat.columns))

regex_e6 = 'B6'
temp_dat = data.iloc[:,data.columns.str.contains(regex_e6)]
temp_dat = temp_dat.iloc[:,~(temp_dat.columns.str.contains('QC'))]
sumzeros = (temp_dat == 0).astype(int).sum(axis=1)
zeroBool_etoh = sumzeros.ge((1-threshold_0)*len(temp_dat.columns))

result = ~(zeroBool_glucose & zeroBool_etoh)
data = data[result == True]

# add on the rows again, and save the output
pd.concat([src_cols.iloc[:,range(0,(src_cols.shape[1]-2))],data],axis = 0)
data.to_excel('data/metabolomics/B2_filtered.xlsx')


# --------- load files, batch1 ---------
src_cols = pd.read_excel('data/metabolomics/Metabolomics_B1.xlsx').iloc[1:4,]
src_data = pd.read_excel('data/metabolomics/Metabolomics_B1.xlsx', sheet_name = 'Full data').iloc[3:,]
src_data.columns = src_data.iloc[0,]
src_data = src_data.iloc[1:,0:(src_data.shape[1]-10)]

# -------- filtering --------
column_name = 'S/N average'
data = SNfilter(src_data,column_name,sn_threshold)

column_name = 'Fill %'
data = Fillfilter(src_data,column_name,fill_threshold)

#filter out if mean is higher than blank?
regex_samples = 'WT|dld3|faa1|gal11|mek1|oca1|pcl1|rme|rts3|rst3|tda1|ygr|cst6'
regex_blanks = 'blank|MeOH'
data = blankfilter(data, regex_blanks, regex_samples)

#filter out zeros in each group
regex_b7 = 'B1G|B2G'
temp_dat = data.iloc[:,data.columns.str.contains(regex_b7)]
temp_dat = temp_dat.iloc[:,~(temp_dat.columns.str.contains('QC'))]
sumzeros = (temp_dat == 0).astype(int).sum(axis=1)
zeroBool_glucose = sumzeros.ge(threshold_0*len(temp_dat.columns))

regex_e6 = 'B1E|B2E'
temp_dat = data.iloc[:,data.columns.str.contains(regex_e6)]
temp_dat = temp_dat.iloc[:,~(temp_dat.columns.str.contains('QC'))]
sumzeros = (temp_dat == 0).astype(int).sum(axis=1)
zeroBool_etoh = sumzeros.ge((1-threshold_0)*len(temp_dat.columns))

result = ~(zeroBool_glucose & zeroBool_etoh)
data = data[result == True]

#now we have filtered data?

# add on the rows again, and save the output
pd.concat([src_cols.iloc[:,range(0,(src_cols.shape[1]-10))],data],axis = 0)
data.to_excel('data/Metabolomics/B1_filtered.xlsx')