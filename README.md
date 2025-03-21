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
> Filename: GeneCountMatrix.csv
>
> Output CSV file structure
> 
> "gene_id","Control_R1","Control_R2","Control_R3","Control_R4","Experimental_R1","Experimental_R2","Experimental_R3","Experimental_R4"
> 
> "gene_id": Ensembl ID
> 
> Columns 2-9: "expected_count" column from respective input files

### **Step 2**

**Filename**: Outlier_Sample_Assessment.R

**Purpose**: Analyze each individual samples and assess whether a sample is an outlier to be removed or not.

**Input**: GeneCountMatrix.csv from step 1

**Output**: PCA plot / Heatmap with batch effect visualization / Indication of outlier samples
