import pandas as pd
import numpy as np
import seaborn as sns
import os
os.chdir ("C:\\Users\\Dell\\Desktop\\Pupil Bio_Tasks\\Task_1")
df = pd.read_csv ('PupilBioTest_PMP_revA.csv')

#PMP's identification with high specificity
df.columns = df.columns.str.replace("`", "")
df['Total_Coverage'] = df[['000', '001', '010', '011', '100', '101', '110', '111']].sum(axis=1)
total_tissue1 = df[df['Tissue'] == 'Tissue1']['Total_Coverage'].sum()
total_tissue2 = df[df['Tissue'] == 'Tissue2']['Total_Coverage'].sum()

from scipy.stats import fisher_exact
def calculate_p_value(row):
    #counts for each tissue
    tissue1_count = row['Total_Coverage'] if row['Tissue'] == 'Tissue1' else 0
    tissue2_count = row['Total_Coverage'] if row['Tissue'] == 'Tissue2' else 0
    
    # Build contingency table
    table = [[tissue1_count, tissue2_count], 
             [total_tissue1 - tissue1_count, total_tissue2 - tissue2_count]]
    
    # Perform Fisher's exact test
    return fisher_exact(table)[1]
    
    #Calculate specificity p-values
df['PMP_Specificity'] = df.apply(calculate_p_value, axis=1)
print("Top 10 Specific PMPs:")
print(df.sort_values('PMP_Specificity').head(10))

df['VRF'] = df[['000', '001', '010', '011', '100', '101', '110', '111']].max(axis=1) / df['Total_Coverage']
vrf_stats = df.groupby('Tissue')['VRF'].mean()
print("Mean VRF for Each Tissue:")
print(vrf_stats)

print(len(df['Total_Coverage']))
print(len(df['PMP_Specificity']))

#Scatterplot to check how sequencing depth can affect specificity
sns.scatterplot(x=df['Total_Coverage'], y=['PMP_Specificity'])
plt.title('Sequencing Depth VS Specificity Confidence')
plt.xlabel('Sequencing Depth')
plt.ylabel('Specificity (p-value)')
plt.yscale('log')
plt.show()

#Estimation of threshol of reads
thr = df[df['Tissue'] == 'Tissue2']['Total_Coverage'].quantile(0.95)
print(f"Threshold of reads required for confident calling of Tissue #2 at a sequencing depth of 1 million reads: {thr}")

#Validation Hypothesis
pmp_specificity = df.groupby('CpG_Coordinates')['PMP_Specificity'].mean()
individual_cpg_specificity = df.groupby(['000', '001', '010', '011', '100', '101', '110', '111'])['Total_Coverage'].mean()

sns.histplot(pmp_specificity, bins=30, label='PMP Specificity', kde=True)
sns.histplot(individual_cpg_specificity, bins=30, label='Individual CpG Specificity', kde=True)
plt.title('Specificity Comparison: PMPs vs Individual CpG Sites')
plt.xlabel('Specificity (p-value)')
plt.ylabel('Frequency')
plt.legend()
plt.show()