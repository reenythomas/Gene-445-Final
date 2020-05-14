import os
os.chdir('/Users/reenythomas/Desktop')

import pandas as pd
import numpy as np

##############
#Objective: To import excel sheet into the work environment
##############

data = pd.read_table('GSE1832_series_matrix-clean.txt', sep='\t', header=0)
data = pd.DataFrame(data) #To change it to a dataframe. 

##############
#Objective: To find the name of the columns 
##############
list(data.columns.values)
    
'''
##############
#Output
##############
!Sample_title
Patient 1 Nonexercise 8 pm
Patient 2 Nonexericse 8 pm
Patient 3 Nonexercise 8 pm
Patient 4 Nonexercise 8 pm
Patient 1 Exercise 8 pm
Pateint 2 Exercise 8 pm
Patient 3 Exercise 8 pm
Patient 4 Exercise 8 pm
Patient 1 Nonexercised 8 am
Patient 2 Nonexercised 8 am
Patient 3 Nonexercised 8 am
Patient 4 Nonexercised 8 am
Patient 2 Exercise 8 am
Patient 3 Exercise 8 am
Patient 4 Exercise 8 am
'''
##############
#Objective: Make the first column the row name
##############    

data = data.set_index('!Sample_title')
data.index.names = [None]

'''
##############
#Output
##############
['Patient 1 Nonexercise 8 pm',
 'Patient 2 Nonexericse 8 pm',
 'Patient 3 Nonexercise 8 pm',
 'Patient 4 Nonexercise 8 pm',
 'Patient 1 Exercise 8 pm',
 'Patient 2 Exercise 8 pm',
 'Patient 3 Exercise 8 pm',
 'Patient 4 Exercise 8 pm',
 'Patient 1 Nonexercised 8 am',
 'Patient 2 Nonexercised 8 am',
 'Patient 3 Nonexercised 8 am',
 'Patient 4 Nonexercised 8 am',
 'Patient 2 Exercise 8 am',
 'Patient 3 Exercise 8 am',
 'Patient 4 Exercise 8 am']
'''

##############
#Objective: Condense the data by log2 transformation
##############    
data = np.log2(data)

##############
#Objective: Calculate F-test
##############  
ftest_list_pm = []
ftest_list_am = []
p_value_list_pm = [] 
p_value_list_am = [] 

import scipy.stats as stats
for i in range(len(data)):
     ftest_pm, p_value_pm = stats.f_oneway(data.iloc[i,0:3], data.iloc[i,4:7]) 
     ftest_am, p_value_am = stats.f_oneway(data.iloc[i,8:11], data.iloc[i,12:14]) 
     ftest_list_pm.append(ftest_pm)
     ftest_list_am.append(ftest_am)
     p_value_list_pm.append(p_value_pm)
     p_value_list_am.append(p_value_am)
     i+=1
      
##############
#Objective: Average for each group: Nonexercise 8 pm, Exercise 8 pm, Nonexercise 8 am, Exercise 8 am
##############  
avg_non8pm_list = []
avg_exe8pm_list = []
avg_non8am_list = []
avg_exe8am_list = []

for i in range(len(data)): 
    avg_non8pm = sum(data.iloc[i,0:3])/len(data.iloc[i,0:3])
    avg_exe8pm = sum(data.iloc[i,4:7])/len(data.iloc[i,4:7])
    avg_non8am = sum(data.iloc[i,8:11])/len(data.iloc[i,8:11])
    avg_exe8am = sum(data.iloc[i,12:14])/len(data.iloc[i,12:14])
    avg_non8pm_list.append(avg_non8pm)
    avg_exe8pm_list.append(avg_exe8pm)
    avg_non8am_list.append(avg_non8am)
    avg_exe8am_list.append(avg_exe8am)
    i+=1

    
##############
#Objective: Calculate log differences 
############## 
diff_8pm_list = []
diff_8am_list = []

for i in range(len(data)): 
    diff_8pm = avg_exe8pm_list[i] - avg_non8pm_list[i]
    diff_8am = avg_exe8am_list[i] - avg_non8am_list[i]
    diff_8pm_list.append(diff_8pm)
    diff_8am_list.append(diff_8am)
    i+=1
    
diff_non_list = []
diff_exe_list = []
    
for i in range(len(data)):
    diff_non = avg_non8pm_list[i] - avg_non8am_list[i]
    diff_exe = avg_exe8pm_list[i] - avg_exe8am_list[i]
    diff_non_list.append(diff_non)
    diff_exe_list.append(diff_exe)
    
#######################
#Objective: Plot Histogram
#######################
import matplotlib.pyplot as plt

plt.hist(diff_8pm_list, bins = 10)
plt.show()  
 
plt.hist(diff_8am_list, bins = 10)
plt.show()  
  
plt.hist(diff_non_list, bins = 10)
plt.show()

plt.hist(diff_exe_list, bins = 10)
plt.show()
 
    
##############
#Objective: Calculate geometric fold changes 
############## 

    
geofold_change_8pm = []
geofold_change_8am = []
    
for i in range(len(data)):
    if diff_8pm_list[i] <0:
        geofold_change = (-1)/(2**(diff_8pm_list[i]))
        geofold_change_8pm.append(geofold_change)
    else:
       geofold_change = (2)**(diff_8pm_list[i])
       geofold_change_8pm.append(geofold_change)
    i += 1  

for i in range(len(data)):
    if diff_8am_list[i] <0:
        geofold_change = (-1)/(2**(diff_8am_list[i]))
        geofold_change_8am.append(geofold_change)
    else:
       geofold_change = (2)**(diff_8am_list[i])
       geofold_change_8am.append(geofold_change)
    i += 1  
    
#############
#Objective: Plot fold changes
############
     
plt.hist(geofold_change_8pm, bins = 5)
plt.show()

plt.hist(geofold_change_8am, bins = 5)
plt.show()
    


