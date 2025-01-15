import pandas as pd
import numpy as np
import os
import matplotlib.pyplot as plt
import seaborn as sns

os.chdir ("C:\\Users\\Dell\\Desktop\\Pupil Bio_Tasks\\Task_1")
df = pd.read_csv ('PupilBioTest_PMP_revA.csv')

df.columns = df.columns.str.replace("`", "")
df['Total_Coverage'] = df[['000', '001', '010', '011', '100', '101', '110', '111']].sum(axis=1)

#Histogram of Total Coverage
sns.histplot(df['Total_Coverage'], kde=True, bins=10)
plt.title('Distribution of Total Coverage')
plt.xlabel('Total Coverage')
plt.ylabel('Frequency')
plt.show()

coverage_categories = ['000', '001', '010', '011', '100', '101', '110', '111']
category_sums = df[coverage_categories].sum()

#Barplot for coverage per category
sns.barplot(x=category_sums.index, y=category_sums.values)
plt.title('Coverage per Category')
plt.xlabel('Category')
plt.ylabel('Total Coverage')
plt.show()

#Boxplot of Total Coverage
sns.boxplot(x=df['Tissue'], y=df['Total_Coverage'])
plt.title('Boxplot of Total Coverage statistics')
plt.ylabel('Total Coverage')
plt.xlabel('Tissue')
plt.show()