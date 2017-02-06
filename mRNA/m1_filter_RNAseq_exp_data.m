%% filter the TCGA RNAseq rsem normalized data by variance and median
%% then input to lmQCMWorkflow_SCC.m for co-expression cluster finding
%% written by Jie Zhang 11.2016
clear
close all

load('KIRC_RNAseq_RSEM_genes_normalized_all.mat');

GeneSymbol_f= GeneSymbol(median(geneExp,2)>1);
exp_f=geneExp((median(geneExp,2)>1),:);

vari = var(exp_f,0,2);
ind2= find(vari>=quantile(vari,0.30));

geneExp_f2=exp_f(ind2,:);
GeneSymbol_f2= GeneSymbol_f(ind2);

save('filtered_geneExp_KIRC.mat','GeneSymbol_f2','geneExp_f2','sampleID');
