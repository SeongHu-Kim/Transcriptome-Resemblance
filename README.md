# Transcriptome-Resemblance

## **Workflow**
### **Step 1**

**Filename**: CountsFiles_to_CountsMatrix.R

**Purpose**: Merge multiple gene counts data into a single CSV file

**Input**: individual raw gene counts .genes.resutls files
>Control_R1.genes.results
>
>Control_R2.genes.results
>
>Control_R3.genes.results
>
>Control_R4.genes.results
>
>Experimental_R1.genes.results
>
>Experimental_R2.genes.results
>
>Experimental_R3.genes.results
>
>Experimental_R4.genes.results

**Output**: Single gene counts matrix text file
> Filename: GeneCountMatrix_Original.csv
>
> Output CSV file structure
> 
>> "gene_id","gene_symbol","Control_R1","Control_R2","Control_R3","Control_R4","Experimental_R1","Experimental_R2","Experimental_R3","Experimental_R4"
>> 
>> "gene_id": Ensembl ID
>>
>> "gene_symbol": Gene symbols respective to Ensembl ID
>>
>> Columns 3-10: "expected_count" column from respective input files

### **Step 2**

**Filename**: Outlier_Sample_Assessment.R

**Purpose**: Analyze each individual samples and assess whether a sample is an outlier to be removed or not.

**Input**: GeneCountMatrix_Original.csv from step 1

**Output**: PCA plot / Heatmap with batch effect visualization / Indication of outlier samples

### **[Optional] Step 2.5**

**Filename**: Remove_Outlier_Samples.R

**Purpose**: Remove outlier samples from GeneCountMatrix_Original.csv

**Input**: GeneCountMatrix_Original.csv from step 1

**Output**: GeneCountMatrix_filtered.csv

### **Step 3**

**Filename**: Differential_Expression_Analysis.R

**Purpose**: Conduct differential expression analysis

**Input**: GeneCountMatrix_Original.csv from step 1 (GeneCountMatrix_filtered.csv from step 2.5 if some samples were dropped)

**Output**: DESeq2 based Differential expression analysis results
> Filename: DESeq2_results.csv, UpDEG.csv, DownDEG.csv, DESeq2_normalized_counts.csv
> 
>> DESeq2_results.csv: Contains all genes
>> 
>> UpDEG.csv: Contains genes with 1.5 fold upregulation and with padj < 0.05
>> 
>> DownDEG.csv: Contains genes with 1.5 fold downregulation and with padj < 0.05
>>
>> DESeq2_normalized_counts.csv: Contains DESeq2 normalized counts for all samples
>> 
> Output CSV file structure (DESeq2_results.csv, UpDEG.csv, DownDEG.csv)
>> 
>> "gene_id","gene_symbol","baseMean","log2FoldChange","lfcSE","stat","pvalue","padj"
>>
>> "gene_id": Ensembl ID
>> 
>> "gene_symbol": Gene symbols respective to Ensembl ID
>>
>> "baseMean": The average normalized count for the gene across all samples (both Control and Experimental)
>>
>> "log2FoldChange": The estimated change in gene expression between the two groups ('Experimental' vs 'Control')
>>
>> "lfcSE": Log Fold Change Standard Error. A smaller lfcSE indicates a more precise estimate of the fold change
>>
>> "stat": The Wald statistic. Larger absolute values of the stat suggest stronger evidence against the null hypothesis
>>
>> "pvalue": The raw p-value obtained from the Wald test. Lower p-values suggest stronger evidence against the null hypothesis
>>
>> "padj": The adjusted p-value (also known as the False Discovery Rate or FDR). Genes with a low padj are considered significantly differentially expressed between the conditions
>>
> Output CSV file structure (DESeq2_normalized_counts.csv)
>
>> Idential to GeneCountMatrix_Original.csv from step 1 (GeneCountMatrix_filtered.csv from step 2.5 if some samples were dropped)

### **[Optional] Step 4A**

**Filename**: Volcano_Plot.R

**Purpose**: Create a Volcano plot based on the DESeq2 results

**Input**: DESeq2_results.csv from step 3

**Output**: volcano_plot.jpg

### **Step 5**

**Filename**: Public_Dataset_to_CountsMatrix.R

**Purpose**: Merge multiple gene counts data into a single CSV file

**Input**: individual raw gene counts .txt files
> Using public available dataset GSE97263

**Output**: Single gene counts matrix text file
> Filename: GeneCountMatrix_GSE97263.csv
>
> Output CSV file structure
> 
>> "gene_id","SLE active_CDR-2","SLE active_CDR-6","SLE active_CDR-20",...
>> 
>> "gene_id": Ensembl ID
>>
>> Columns 2-45: raw gene counts column from respective input files
>
### **Step 6**

**Filename**: Combine_CountsMatrix.R

**Purpose**: Merge personal and public gene counts matrices

**Input**: Personal and public gene counts matrices
> Personal gene counts matrix: GeneCountMatrix_Original.csv from step 1 (GeneCountMatrix_filtered.csv from step 2.5 if some samples were dropped)
> 
> Public gene counts matrix: GeneCountMatrix_GSE97263.csv from step 5

**Output**: Single combined gene counts matrix text file
> Filename: GeneCountMatrix_Combined.csv
>
> Output CSV file structure
> 
>> "gene_id","gene_symbol","Control_R1","Control_R3",...
>> 
>> "gene_id": Ensembl ID
>>
>> "gene_symbol": Gene symbols respective to Ensembl ID
>>
>> Columns 3-52: raw gene counts column from respective input files

