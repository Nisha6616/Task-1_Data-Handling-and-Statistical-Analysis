import pandas as pd
import numpy as np
import os
import statistics
#Load the dataset
os.chdir ("C:\\Users\\Dell\\Desktop\\Pupil Bio_Tasks\\Task_1\\Sub_A")
df = pd.read_csv ('PupilBioTest_PMP_revA.csv')

# Remove backticks from column names
df.columns = df.columns.str.replace("`", "")
df1 = ['000', '001', '010', '011', '100', '101', '110', '111']
# separate DataFrames for 'cfDNA' and 'Islet' tissues
cfdna_df = df[df['Tissue'] == 'cfDNA']
islet_df = df[df['Tissue'] == 'Islet']
cfdna_output_file = "cfdna.csv"
islet_output_file = "islet.csv"
cfdna_df.to_csv(cfdna_output_file, index=False)
islet_df.to_csv(islet_output_file, index=False)

#P ROCESSING CFDNA DATA
new1 = pd.read_csv ('cfdna.csv')
new1 = df.assign(Total=new1['000']+new1['001']+new1['010']+new1['011']+new1['100']+new1['101']+new1['110']+new1['111'])
print (new1.head(10))

median = new1['Total'].median()
print ("Median is ",median)
mean = new1['Total'].mean()
print ("Mean is ",mean)
length = len(new1['Total'])
print ("Length is ",length)

# coefficient of variation (CV) calculation
squared_diffs = [(x - mean) ** 2 for x in new1['Total']]
cv_coverage = (std_dev / mean) * 100
print (cv_coverage)

# PROCESSING ISLET DATA
islet_df.to_csv(islet_output_file, index=False)
new2 = pd.read_csv ('islet.csv')
new2 = df.assign(Total=new2['000']+new2['001']+new2['010']+new2['011']+new2['100']+new2['101']+new2['110']+new2['111'])
median = new2['Total'].median()
print ("Median is ",median)
mean = new2['Total'].mean()
print ("Mean is ",mean)
length = len(new2['Total'])
print ("Length is ",length)

# coefficient of variation (CV) calculation
std_dev = np.std(new2['Total'])
cv_coverage = (std_dev / mean) * 100
print (cv_coverage)