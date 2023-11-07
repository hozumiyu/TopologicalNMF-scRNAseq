# TopologicalNMF-scRNAseq
## This repository can be used to produce the results for the paper entitled 'Analyzing Single Cell RNA Sequencing with Topological Nonnegative Matrix Factorization'

## Content:
* main_nmf_benchmark.py
* * This produces the result for NMF, rNMF, GNMF, rGNMF benchmark results
  * run python main_nmf_benchmark.py dataset
  * dataset = 'GSE67835', 'GSE64016', 'GSE75748time', 'GSE82187', 'GSE84133human1', 'GSE84133human2', 'GSE84133human3', 'GSE84133human4', 'GSE84133mouse1', 'GSE84133mouse2', 'GSE57249', 'GSE94820'
* main_tnmf.py
* * This produces the results form topological nmf
  * python main_tnmf.py dataset method n_neighbors weights l
  * method = 'TNMF', 'kTNMF', 'rTNMF', 'krTNMF'
  * For the paper, n_neighbors = 8 and l = 1 was used.
  * weights must be delimitted by comma, for example 1,0,0,1,0,1
  * example: python main_tnmf.py GSE67835 TNMF 8 1,0,0,1,0,1 1
