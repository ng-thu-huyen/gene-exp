# gene-exp
The dataset using in this project is from the open access database on GTEx Portal.  
- Gene expression: https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SampleAttributesDD.xlsx  
- Sample attribute: https://storage.googleapis.com/gtex_analysis_v8/rna_seq_data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.gct.gz  
- Subject phenotypes: https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt  

## convert gene expression .gct file to DataFrame
- using the convertgct file
```
import pandas as pd
with open(file_path, 'r') as f:
  lines = f.readlines()
cols = lines[2].split()
raw_df = {col: [] for col in cols}
for idx, row in enumerate(lines[3:]):
  cells = row.split()
  if len(cells) < 17384:
    continue
  for i, col in enumerate(cols):
    raw_df[col].append(cells[i])
df = pd.DataFrame(raw_df)
df.to_csv(file_path_to_save)
```

## filter only whole blood sample IDs
- find all whole blood sampld IDs in Sample attribute file where column 'SMTSD' = 'Whole Blood'

```
import pandas as pd
df = pd.read_csv('https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SampleAttributesDS.txt', delimiter = "\t")
blood_df = df.loc[df['SMTSD'] == "Whole Blood"] 
whole_samp = blood_df['SAMPID'] #list that contains whole blood sample ids
```

- filter whole blood sample ids in the gene expression file
```
exp_samp = gene.columns
sampleinter = ['Description']
for sampleid1 in whole_samp:
  for sampleid2 in exp_samp:
    if sampleid1 == sampleid2:
      sampleinter.append(sampleid1)
gene_1 = gene[sampleinter]
gene_1
```
## filter geneids that already have coefficients (based on the formula on the [paper](https://www.nature.com/articles/ncomms9570#MOESM435))
```
import pandas as pd
#read geneid-coefficient file
co = pd.read_csv(file_path)
co = co.loc[:, ~co.columns.str.contains('^Unnamed')]
#there are duplicate geneids with different coefficient, only keep the first row
co = co.drop_duplicates(subset='geneid', keep="first")
#dict with key is geneid and value is its coefficient
codict = {key:value for key, value in zip(co['geneid'], co['factor'])}
sub_gene = gene_1.loc[gene_1['Description'].isin(co['geneid'])]
sub_gene
```
## calculate predicted age (based on the formula on the [paper](https://www.nature.com/articles/ncomms9570#MOESM435))
```
import numpy
dict_age = {}
for col in subgene.columns[1:-1]:
  a = subgene[col].to_numpy()
  b = subgene['factor'].to_numpy()
  age = numpy.dot(subgene[col].to_numpy(), subgene['factor'].to_numpy())
  dict_age[col[0:-14]] = age
age = pd.DataFrame(dict_age.items(), columns = ['SUBJID', 'PredictedAge'])
#open Subject phenotypes file (for real age of subject)
real_age = pd.read_csv('https://storage.googleapis.com/gtex_analysis_v8/annotations/GTEx_Analysis_v8_Annotations_SubjectPhenotypesDS.txt', delimiter = "\t")
agedict = {key:value for key, value in zip(real_age['SUBJID'], real_age['AGE'])}
age['RealAge'] = age['SUBJID'].apply(lambda g: agedict[g])
age
```
## draw boxplot
```
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt 

age.boxplot(by ='RealAge', column = 'PredictedAge', grid = False)
plt.savefig('Boxplot(755 samples)')
```
## initial findings
- The paper uses 11779 gene ids and 7074 whole blood samples 
- Using GTEx Portal v8, we have 56203 gene ids and 755 whole blood samples 
- The current stats is analyzed based on 11188 intersected gene ids

![Boxplot using v8 GTEx Portal](https://github.com/ng-thu-huyen/gene-exp/blob/main/Boxplot(755%20samples).png))



