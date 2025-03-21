# Transcriptome-Resemblance

## **Workflow**
### **Step 1**

Filename: CountsFiles_to_CountsMatrix.R

Purpose: Merge multiple gene counts data into a single CSV file

Input: individual raw gene counts .genes.resutls files

Output: single GeneCountMatrix.csv file
> Output CSV file structure
> 
> "gene_id","Control_R1","Control_R2","Control_R3","Control_R4","Experimental_R1","Experimental_R2","Experimental_R3","Experimental_R4"
> 
> "gene_id": Ensembl ID
> 
> Columns 2-9: "expected_count" column from respective input files

