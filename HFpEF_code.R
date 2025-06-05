#####HFpEF Survival Analysis#####
#################################
#Version: 3.2
#Date: June 05, 2025
#Author: Suraj Kannan (working with Dr. Virginia Hahn)

library(Matrix)
library(DESeq2)
library(pheatmap)
library(RColorBrewer)
library(randomForestSRC)
library(stringr)
library(survival)
library(survminer)
library(survMisc)
library(data.table)
library(purrr)
library(dplyr)
library(clipr)
library(ggrepel)
library(sceasy)
library(reticulate)
library(Seurat)

setwd("~/Documents/Research/HFpEF/")
fig_factor = 1

data_summary <- function(data, varname, groupnames){
  require(plyr)
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}

###Data Loading###
#The first workspace contains all of the RNA-seq data; we will do some clean-up here to get it into a consistent format
load("vsh_getcounts.RData")
samplesInfo = samplesInfo[!samplesInfo$SequenceRunId %in% excludeInfo$SequenceRunId, ] #there are some samples excluded because of clinical findings
samplesTable$subjectid = samplesInfo$`Subject ID`
samplesTable$ethnicity = samplesInfo$`Subject Ethnicity`
#I'm going to rename some objects to fit my standard conventions
data = counts
pheno = samplesTable
#Now we need to appropriately format the data and pheno tables
rownames(pheno) = paste("s", pheno$samplenames, sep = "")
colnames(data) = paste("s", colnames(data), sep = "")
data = as.data.frame(data)
#As a last piece of formatting, we need to convert the ensembl gene names to symbols just for ease of working. There are multiple ways to do this but I will take the laziest approach and use the code I've already written for this purpose in another project.
load("~/Documents/Research/BigRef/FinalWorkingFiles/clean_nodatasets_060720.RData")
data = rename_genes(data, species = "human") #Note that doing this cuts down the number of genes pretty significantly (from 50k -> 30k), but keeps about 99.4% of the counts. So it should be fine!
data = as.data.frame(as.matrix(data))
rm(list=setdiff(ls(), c("data", "pheno", "fig_factor", "data_summary"))) #The only objects we need are the data and pheno objects so we can get rid of everything else
save_data = data #turns out that there's an object in the hospitalization data called "data" so this is a clunky intermediate

#Now we will load in and format the hospitalization data
load("2020rnaseq_hfpefsurvival.RData")
pheno_hfpef = pheno[pheno$disease == "HFpEF", ]
#Our pheno table, any_hosp, cv_hosp, and hf_hosp are all labeled in different orders, so let's fix that
any_hosp = as.data.frame(any_hosp)
cv_hosp = as.data.frame(cv_hosp)
hf_hosp = as.data.frame(hf_hosp)
rownames(any_hosp) = paste("s", any_hosp$SequenceRunId, sep ="")
rownames(cv_hosp) = paste("s", cv_hosp$SequenceRunId, sep ="")
rownames(hf_hosp) = paste("s", hf_hosp$SequenceRunId, sep ="")
any_hosp = any_hosp[rownames(pheno_hfpef), ]
cv_hosp = cv_hosp[rownames(pheno_hfpef), ]
hf_hosp = hf_hosp[rownames(pheno_hfpef), ]

#In the earlier version of this code, we updated the HFpEF pheno table with information from the other tables to note hospitalizations/events. However, we now have an updated table from Virginia that has HF hospitalizations + any death data for our patients, which is actually what we want to be focused on. There are some discrepancies between this recent table and the data in the other loaded tables, but we are going to use the most recent version from Virginia.
virginia_temp = read.csv("survival_nppb.csv", as.is = TRUE, header = TRUE)
rownames(virginia_temp) = paste(virginia_temp$Biopsy.ID, "p", sep = "")
pheno_hfpef$time_to_event = virginia_temp[pheno_hfpef$subjectid, ]$futime_3yhfdeath
pheno_hfpef$event_stat = virginia_temp[pheno_hfpef$subjectid, ]$fustat_3yhfdeath
pheno_hfpef$year = factor(trunc(pheno_hfpef$time_to_event/366) + 1)
pheno_hfpef$category = paste("Year", pheno_hfpef$year)
pheno_hfpef[pheno_hfpef$year == 1 & pheno_hfpef$event_stat == 0, ]$category = "Year 1 Lost"
pheno_hfpef[pheno_hfpef$year == 3, ]$category = "No Event"
pheno_hfpef$category = factor(pheno_hfpef$category, levels = c("Year 1", "Year 1 Lost", "Year 2", "No Event"))
#As a note - in the original analysis, there were 6 patients excluded. Three were excluded because they were RNA-seq outliers for some reason; however, I don't see a reason to exclude them from my analysis. The other three were initially excluded due to early censoring; however, this doesn't appear to be the case in the new updated survival data. As a result, I will be keeping all of the patients until further assessment. 
#Lastly, let's get rid of all the excess objects
data = save_data
rm(list=setdiff(ls(), c("data", "pheno", "any_hosp", "cv_hosp", "hf_hosp", "pheno_hfpef", "fig_factor", "data_summary")))

###HFpEF and HFrEF differential gene expression analysis###
#This largely repeats analysis from the original manuscript; it's not necessary but I'm keeping on hand just in case
dds = DESeqDataSetFromMatrix(data[rowMeans(data) >= 50, pheno$tissue == "RVS"], pheno[pheno$tissue == "RVS", ], design = ~sex + disease) #Select only RVS tisuse and set up DESeq object; I used a mean count of 50 per gene as a minimum, similar to the original paper
dds = DESeq(dds)
res = results(dds)
rld = vst(dds, blind = FALSE)

#Extract out the three comparisons - HFpEF vs control, HFrEF vs control, and HFpEF vs HFrEF
hfpef = as.data.frame(results(dds, contrast = c("disease", "HFpEF", "Normal")))
hfpef = hfpef[order(hfpef$padj), ]
hfref = as.data.frame(results(dds, contrast = c("disease", "HFrEF", "Normal")))
hfref = hfref[order(hfref$padj), ]
hfcomp = as.data.frame(results(dds, contrast = c("disease", "HFpEF", "HFrEF")))
hfcomp = hfcomp[order(hfcomp$padj), ]
hfpef_genes = rownames(hfpef)[hfpef$padj < 0.05]
hfref_genes = rownames(hfref)[hfref$padj < 0.05]
hfcomp_genes = rownames(hfcomp)[hfcomp$padj < 0.05]

###Differential gene expression analysis for survival groups###
#Just as preliminary analysis, I start by splitting patients into "year of relapse" (calculated from time to event - note that patients with "year 3" as their year of relapse never had a HF/all-cause mortality event). For the sake of this analysis, I eliminated the three patients who were lost to follow-up (e.g. time to event within first year but no actual event recorded) - we can of course include those patients later in the random survival forest analysis.

dds_survival = DESeqDataSetFromMatrix(data[rowMeans(data) >= 50, rownames(pheno_hfpef)[!pheno_hfpef$category == "Year 1 Lost"]], pheno_hfpef[rownames(pheno_hfpef)[!pheno_hfpef$category == "Year 1 Lost"], ], design = ~sex + year) #Select only HFpEF samples (again discarding the ones with limited followup) and set up DESeq object; I used a mean count of 50 per gene as a minimum, similar to the original paper
dds_survival = DESeq(dds_survival)
res_survival = results(dds_survival)
survival = as.data.frame(results(dds_survival, contrast = c("year", "1", "3")))
survival = survival[order(survival$padj), ]

#Now, it becomes clear that this method only yields a small number of differentially expressed genes. My suspicion is that this is because the data is somewhat underpowered. An interesting and easy way to further check this is to instead divide across whether the samples are from year 1 or not - this adds the year 2 samples together with the year 3 samples (code for this is below, commented out). You would think that this should change the results too much given there are only five year 2 samples, but the number of differentially expressed genes *does* increase from 15 to 33. Also worth noting that genes like NPPB have an FDR hovering just above 0.05, which I think is pretty indicative of underpowering.

#pheno_hfpef$year1 = (pheno_hfpef$year == 1)
#dds_survival = DESeqDataSetFromMatrix(data[rowMeans(data) >= 50, rownames(pheno_hfpef)[!(pheno_hfpef$year == 1 & pheno_hfpef$event_stat == 0)]], pheno_hfpef[rownames(pheno_hfpef)[!(pheno_hfpef$year == 1 & pheno_hfpef$event_stat == 0)], ], design = ~sex + year1) #Select only HFpEF samples (again discarding the ones with limited followup) and set up DESeq object; I used a mean count of 50 per gene as a minimum, similar to the original paper
#dds_survival = DESeq(dds_survival)
#res_survival = results(dds_survival)
#survival = as.data.frame(results(dds_survival))
#survival = survival[order(survival$padj), ]

#By the way, I think it is worth noting that sex differences really make a huge difference; removing sex as a confounder above increases the number of differentially expressed genes pretty significantly. Thankfully, our random forest approach should be able to handle non-linear effects like that, but it's very much worth keeping in mind.

#Now, let's make a PCA plot of the HFpEF samples.
pca_hfpef = prcomp(assay(rld)[, rownames(pheno_hfpef)])
pca_hfpef_plot = as.data.frame(pca_hfpef$rotation)
pca_hfpef_plot$category = pheno_hfpef$category
ggplot(pca_hfpef_plot, aes(x = PC1, y = PC2, color = category)) + geom_point(size = 4/fig_factor) + theme_linedraw() + theme(axis.title = element_text(size = 16/fig_factor), axis.text = element_text(size = 12/fig_factor), legend.title = element_text(size = 12/fig_factor), legend.text = element_text(size = 10/fig_factor)) + guides(color = guide_legend(title = "Year of Event")) + scale_color_manual(values = c("indianred1", "royalblue", "orange", "gold")) + ylab("PC2") + xlab("PC1")

###HFpEF Random Survival Forest Optimization Steps###
#For survival analysis, I'm using the random survival forest package as well as several associated strategies from Hemant Ishwaran (the most relevant papers are his papers on high-dimensional variable selection for survival). I am identifying the critical HFpEF survival signature by using variable hunting method, which combines stepwise regularization with minimal depth thresholding (the latter being an alternative to VIMP).

#As with any analysis with so many parameters that I can affect the results, I started with a pretty extensive optimization. I'll include some of the important optimization steps here. As a standard upfront starting point, I often sub-selected genes based on variation, with the assumption that genes with a higher standard coefficient of variation across HFpEF samples were more likely to play a role in driving survival differences. Of course, this isn't necessarily true - after all, we already saw that PC1 was pretty heavily dominated by technical genes above! But those genes don't really predict survival and so it works out fine.

#I started by making a data table of variance-stabilized HFpEF data, inluding only genes with sufficient expression. For want of a better name, this ended up being called "clean_data"
clean_data = vst(as.matrix(data[rowMeans(data) > 50, rownames(pheno_hfpef)]))

#This snippet calculates the coefficient of variation for all genes, and can be used to subselect the desired number of genes
testVars = rowSds(as.matrix(clean_data))/rowMeans(as.matrix(clean_data))
rownames(clean_data) = str_replace_all(rownames(clean_data), "-", ".") #Unfortunately, we have to do this now to avoid some issues down the line with gene names
names(testVars) = rownames(clean_data)
testVars = testVars[!is.na(testVars)]
n = 4000 #Set this to whatever is desired; I selected 4000 on a whim here but will show later that this is a reasonable choice
var_genes = names(sort(testVars, decreasing = TRUE))[1:n]

#Now, above, we did differential testing using DESeq2. As an alternative method, we can fit a Cox proportional hazards model to our survival data. Ishwaran et al. describe in their manuscripts why this is not necessarily an optimal approach, but we will use it below for some of forest applications (for weighting selection of genes). So I will include the code to calculate this here (this code is included here because it calculates off of the normalized clean data).
gene_coxph = numeric()
for(gene in rownames(clean_data)){
  print(gene)
  form = as.formula(paste("Surv(time = pheno_hfpef$time_to_event, event = pheno_hfpef$event_stat) ~ ", gene, sep = ""))
  fit.coxph = coxph(form, data = as.data.frame(t(clean_data)))
  gene_coxph[gene] = summary(fit.coxph)$coefficients[5]
}
rm(form, fit.coxph)

#So now we need some generic code to optimize the forests. The following code can be readily adapted to test whatever variables need to be tested. I include some commentary below.

gene_weights = 1/unlist(lapply(gene_coxph[var_genes], function(x) {max(x, n^(-4/5))})) #Selecting our weights; this can be adjusted to whatever method you choose, but as an initial test I'm showing it for the coxph p-values as weighting
forest_data = as.data.frame(t(clean_data[var_genes, ])) #Creating the data object to be used for the random forest
forest_data$time_to_event = pheno_hfpef$time_to_event
forest_data$event_stat = pheno_hfpef$event_stat
#test_var = tune(Surv(time_to_event, event_stat)~., forest_data)
forest = var.select(object = rfsrc(Surv(time_to_event, event_stat)~., forest_data, importance = TRUE), data = forest_data, method = "vh", mvars = 500, conservative = "medium", xvar.wt = gene_weights, nrep = 100, nstep = 5, ntree = 1000, nodesize = 1, mtry = n^(4/5), nsplit = 10)
var_check = as.data.frame(forest$varselect) #Using variable hunting
var_check$topvar = rownames(var_check) %in% forest$topvars
ggplot(var_check, aes(x = seq(1, length(rownames(var_check))), y = depth, color = topvar)) + geom_point(size = 1/fig_factor) + theme_linedraw() + theme(axis.title = element_text(size = 16/fig_factor), axis.text = element_text(size = 12/fig_factor), legend.title = element_text(size = 12/fig_factor), legend.text = element_text(size = 10/fig_factor)) + guides(color = guide_legend(title = "Top Variable?")) + scale_color_manual(values = c("indianred1", "royalblue")) + ylab("Minimal Depth") + xlab("Ranked Gene Number") #Plotting out the selected top variables

#Optimizing number of trees per forest - I did some pilot tests for this and found that honestly, it doesn't matter as long as you have enough trees. You can adjust this by adjusting ntree on line 162 in the code above. I ended up using 1000, which I think was plenty.
#Optimizing the number of iterations of variable hunting - this has more impact, because with p genes >>> n samples, you could simply not get enough data to assess the variables (they don't get sampled in enough runs) to determine an accurate importance. You can adjust this by adjusting nrep on line 158 in the code above.
#Adjusting tree depth - this does have some impact, because the deeper your trees are, the more range of minimal depth you have to select out the best genes. Ishwaran et al. in general recommend trying to grow the largest trees possible; however, they also have some code to "optimize" the tree depth by creating sample trees with a range of nodesizes and mtrys (commented out on line 161 - the resultant test_var$optimal will have the optimized values). However, my experience was that these select more false positive genes than just manually setting the nodesize to 1 for the deepest trees. mtry doesn't have as much effect, but Ishwaran et al. recommends using a large nodesize, you I just picked n^(4/5) as they do in several of their papers.

#Optimizing the number of input genes n and the method for weighting - so now, we optimize the number of genes as input to the forest, and how we weight. Ishwaran et al. do recommend weighting so I tried four options - uniform (random) weights, weights based on Cox PH modeling (as above), weights based on DESeq2 (I used patients with events in 1 year vs patients with no events), and VIMP (built into random forests). The optimization code is below. Please note - THIS CODE TAKES FOREVER (multiple days). Rather than run it again, I've made a dataspace object that I just load in with all of the relevant objects - the actual code is commented out.

gene_number = c(500, 750, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 5000, 6000, 7000, 8000, 10000)
# list_uniform = list()
# list_coxph = list()
# list_deseq2 = list()
# list_vimp = list()
# for (n in gene_number){
#   var_genes = names(sort(testVars, decreasing = TRUE))[1:n]
#   gene_weights_uniform = rep(1, 1, n)
#   names(gene_weights_uniform) = var_genes
#   gene_weights_coxph = 1/unlist(lapply(gene_coxph[var_genes], function(x) {max(x, n^(-4/5))}))
#   names(gene_weights_coxph) = var_genes
#   gene_weights_deseq2 = 1/unlist(lapply(survival[var_genes, ]$padj, function(x) {max(x, n^(-4/5))}))
#   forest_data_var = as.data.frame(t(clean_data[var_genes, ]))
#   forest_data_var$time_to_event = pheno_hfpef$time_to_event
#   forest_data_var$event_stat = pheno_hfpef$event_stat
#   print(n)
#   list_uniform[[n]] = var.select(object = rfsrc(Surv(time_to_event, event_stat)~., forest_data_var, importance = TRUE), data = forest_data_var, method = "vh", mvars = 500, conservative = "medium", xvar.wt = gene_weights_uniform, nrep = 500, nstep = 5, ntree = 1000, nodesize = 1, mtry = n^(4/5), nsplit = 10, verbose = FALSE)
#   print("uniform")
#   list_coxph[[n]] = var.select(object = rfsrc(Surv(time_to_event, event_stat)~., forest_data_var, importance = TRUE), data = forest_data_var, method = "vh", mvars = 500, conservative = "medium", xvar.wt = gene_weights_coxph, nrep = 500, nstep = 5, ntree = 1000, nodesize = 1, mtry = n^(4/5), nsplit = 10, verbose = FALSE)
#   print("coxph")
#   list_deseq2[[n]] = var.select(object = rfsrc(Surv(time_to_event, event_stat)~., forest_data_var, importance = TRUE), data = forest_data_var, method = "vh", mvars = 500, conservative = "medium", xvar.wt = gene_weights_deseq2, nrep = 500, nstep = 5, ntree = 1000, nodesize = 1, mtry = n^(4/5), nsplit = 10, verbose = FALSE)
#   print("deseq2")
#   list_vimp[[n]] = var.select(object = rfsrc(Surv(time_to_event, event_stat)~., forest_data_var, importance = TRUE), data = forest_data_var, method = "vh", mvars = 500, conservative = "medium", nrep = 500, nstep = 5, ntree = 1000, nodesize = 1, mtry = n^(4/5), nsplit = 10, verbose = FALSE)
#   print("vimp")
# }
# 
# errors_uniform = numeric(14)
# for(n in seq(1,14)){
#   errors_uniform[n] = get.cindex(list_uniform[[gene_number[n]]][["rfsrc.refit.obj"]][["yvar"]][,1], list_uniform[[gene_number[n]]][["rfsrc.refit.obj"]][["yvar"]][,2], list_uniform[[gene_number[n]]][["rfsrc.refit.obj"]][["predicted.oob"]])
# }
# errors_coxph = numeric(14)
# for(n in seq(1,14)){
#   errors_coxph[n] = get.cindex(list_coxph[[gene_number[n]]][["rfsrc.refit.obj"]][["yvar"]][,1], list_coxph[[gene_number[n]]][["rfsrc.refit.obj"]][["yvar"]][,2], list_coxph[[gene_number[n]]][["rfsrc.refit.obj"]][["predicted.oob"]])
# }
# errors_deseq2 = numeric(14)
# for(n in seq(1,14)){
#   errors_deseq2[n] = get.cindex(list_deseq2[[gene_number[n]]][["rfsrc.refit.obj"]][["yvar"]][,1], list_deseq2[[gene_number[n]]][["rfsrc.refit.obj"]][["yvar"]][,2], list_deseq2[[gene_number[n]]][["rfsrc.refit.obj"]][["predicted.oob"]])
# }
# errors_vimp = numeric(14)
# for(n in seq(1,14)){
#   errors_vimp[n] = get.cindex(list_vimp[[gene_number[n]]][["rfsrc.refit.obj"]][["yvar"]][,1], list_vimp[[gene_number[n]]][["rfsrc.refit.obj"]][["yvar"]][,2], list_vimp[[gene_number[n]]][["rfsrc.refit.obj"]][["predicted.oob"]])
# }
# 
# length_uniform = numeric(14)
# for(n in seq(1,14)){
#   length_uniform[n] = length(list_uniform[[gene_number[n]]][["topvars"]])
# }
# length_coxph = numeric(14)
# for(n in seq(1,14)){
#   length_coxph[n] = length(list_coxph[[gene_number[n]]][["topvars"]])
# }
# length_deseq2 = numeric(14)
# for(n in seq(1,14)){
#   length_deseq2[n] = length(list_deseq2[[gene_number[n]]][["topvars"]])
# }
# length_vimp = numeric(14)
# for(n in seq(1,14)){
#   length_vimp[n] = length(list_vimp[[gene_number[n]]][["topvars"]])
# }
# 
# cor_uniform = numeric(14)
# for(n in seq(1,14)){
#   cor_uniform[n] = cor(list_uniform[[gene_number[n]]][["rfsrc.refit.obj"]][["predicted"]][!(forest_data_var$time_to_event < 265 & forest_data_var$event_stat == 0)], forest_data_var$time_to_event[!(forest_data_var$time_to_event < 265 & forest_data_var$event_stat == 0)])
# }
# cor_coxph = numeric(14)
# for(n in seq(1,14)){
#   cor_coxph[n] = cor(list_coxph[[gene_number[n]]][["rfsrc.refit.obj"]][["predicted"]][!(forest_data_var$time_to_event < 265 & forest_data_var$event_stat == 0)], forest_data_var$time_to_event[!(forest_data_var$time_to_event < 265 & forest_data_var$event_stat == 0)])
# }
# cor_deseq2 = numeric(14)
# for(n in seq(1,14)){
#   cor_deseq2[n] = cor(list_deseq2[[gene_number[n]]][["rfsrc.refit.obj"]][["predicted"]][!(forest_data_var$time_to_event < 265 & forest_data_var$event_stat == 0)], forest_data_var$time_to_event[!(forest_data_var$time_to_event < 265 & forest_data_var$event_stat == 0)])
# }
# cor_vimp = numeric(14)
# for(n in seq(1,14)){
#   cor_vimp[n] = cor(list_vimp[[gene_number[n]]][["rfsrc.refit.obj"]][["predicted"]][!(forest_data_var$time_to_event < 265 & forest_data_var$event_stat == 0)], forest_data_var$time_to_event[!(forest_data_var$time_to_event < 265 & forest_data_var$event_stat == 0)])
# }
# 
# errors = data.frame(method = c(rep("uniform", 14), rep("coxph", 14), rep("deseq2", 14), rep("vimp", 14)), genes = rep(gene_number, 4), error = c(errors_uniform, errors_coxph, errors_deseq2, errors_vimp), length = c(length_uniform, length_coxph, length_deseq2, length_vimp), cor = c(cor_uniform, cor_coxph, cor_deseq2, cor_vimp))

#The most important objects from the above code as the "errors" dataframe, and the respective "list" objects. To save a ton of headache, I will load them in here.
load("hfpef_optimization.RData")
ggplot(errors, aes(x = genes, y = error, color = method)) + geom_point(size = 4/fig_factor) + geom_line() + theme_linedraw() + theme(axis.title = element_text(size = 16/fig_factor), axis.text = element_text(size = 12/fig_factor), legend.title = element_text(size = 12/fig_factor), legend.text = element_text(size = 10/fig_factor)) + guides(color = guide_legend(title = "Weighting Method")) + scale_color_manual(values = c("indianred1", "royalblue", "forestgreen", "gold"), labels = c("Cox Proportional Hazard", "DESeq2", "Uniform/Random", "VIMP")) + ylab("Forest Error") + xlab("Number of Input Genes") #Plotting

#You will note that in the above code, I frequently used the "medium" conservativeness. But there is also an option to go more liberal. Using the optimized n = 4000, I calculated the liberal lists for each of the above methods. Again, this is loaded in through the workspace above, so you don't need to re-run it.

# n = 4000
# var_genes = names(sort(testVars, decreasing = TRUE))[1:n]
# gene_weights_uniform = rep(1, 1, n)
# names(gene_weights_uniform) = var_genes
# gene_weights_coxph = 1/unlist(lapply(gene_coxph[var_genes], function(x) {max(x, n^(-4/5))}))
# names(gene_weights_coxph) = var_genes
# gene_weights_deseq2 = 1/unlist(lapply(survival[var_genes, ]$padj, function(x) {max(x, n^(-4/5))}))
# forest_data_var = as.data.frame(t(clean_data[var_genes, ]))
# forest_data_var$time_to_event = pheno_hfpef$time_to_event
# forest_data_var$event_stat = pheno_hfpef$event_stat
# liberal = list()
# print("Start")
# liberal[["uniform"]] = var.select(object = rfsrc(Surv(time_to_event, event_stat)~., forest_data_var, importance = TRUE), data = forest_data_var, method = "vh", mvars = 500, conservative = "low", xvar.wt = gene_weights_uniform, nrep = 500, nstep = 5, ntree = 1000, nodesize = 1, mtry = n^(4/5), nsplit = 10, verbose = FALSE)
# print("Uniform")
# liberal[["coxph"]] = var.select(object = rfsrc(Surv(time_to_event, event_stat)~., forest_data_var, importance = TRUE), data = forest_data_var, method = "vh", mvars = 500, conservative = "low", xvar.wt = gene_weights_coxph, nrep = 500, nstep = 5, ntree = 1000, nodesize = 1, mtry = n^(4/5), nsplit = 10, verbose = FALSE)
# print("Coxph")
# liberal[["deseq2"]] = var.select(object = rfsrc(Surv(time_to_event, event_stat)~., forest_data_var, importance = TRUE), data = forest_data_var, method = "vh", mvars = 500, conservative = "low", xvar.wt = gene_weights_deseq2, nrep = 500, nstep = 5, ntree = 1000, nodesize = 1, mtry = n^(4/5), nsplit = 10, verbose = FALSE)
# print("DESeq2")
# liberal[["vimp"]] = var.select(object = rfsrc(Surv(time_to_event, event_stat)~., forest_data_var, importance = TRUE), data = forest_data_var, method = "vh", mvars = 500, conservative = "low", nrep = 500, nstep = 5, ntree = 1000, nodesize = 1, mtry = n^(4/5), nsplit = 10, verbose = FALSE)
# print("VIMP")

###Finalizing the optimization###
#What we have done above was to identify optimal parameters for the random survival forest (particularly, the number of trees, the number of iterations, and the tree depth), and then run the random survival forest using those parameters at multiple different numbers of input genes using different weighting strategies. Each of those approaches yields an identified survival gene list based on the protocol used by the variable hunting method (again, see the manuscripts by Hemant Ishwaran). What we need to do now is come up with an optimized final gene list from all of those sample tests. As a first step, what we will do is create survival forests based on the identified surivival gene lists for each combination done above (note - we don't need to re-run variable hunting here since we have already hunted the variables! Now we are seeing which identified list is best).

#The premise that I will come to (once all of the code below is run) is that no particular weighting method is necessarily best, and that it may instead be best to pool information from all of the weightings. You can imagine two ways to do this - by taking all the genes that appear across the methods (a union) or by taking genes identified through multiple weighting strategies (intersection). I do both below (with the caveat that the intersection list is generated as genes that are identified by at least three of the four weighting methods for any given input gene number).

first_round_errors = errors #I accidentally used the name `errors` for two different dataframes, so we'll correct this one

#We create the union list and the intersection list at each of the tested input gene numbers
union_list = list()
intersection_list = list()
for(n in gene_number){
  print(n)
  combined = c(list_coxph[[n]][["topvars"]], list_deseq2[[n]][["topvars"]], list_uniform[[n]][["topvars"]], list_vimp[[n]][["topvars"]])
  tabled = table(combined)
  union_list[[n]] = unique(combined)
  intersection_list[[n]] = names(tabled)[tabled >= 3]
  rm(combined, tabled)
}

#What we do here is run the random forest for the identified gene list for each weighting at each gene number, as well as for the union and intersection lists at that gene number. Again - we don't need to rerun variable hunting here, so this is much faster. I ran 5 trials for each method so that we can get a bit of a range.
errors = data.frame()
for(n in gene_number){
  print(n)
  var_genes = names(sort(testVars, decreasing = TRUE))[1:n]
  forest_data = as.data.frame(t(clean_data[var_genes, ])) #Creating the data object to be used for the random forest
  forest_data$time_to_event = pheno_hfpef$time_to_event
  forest_data$event_stat = pheno_hfpef$event_stat
  
  for(trial in seq(1, 5)){
    print(paste("Cox PH", trial))
    obj = rfsrc(Surv(time_to_event, event_stat)~., forest_data[, c(list_coxph[[n]][["topvars"]], "time_to_event", "event_stat")], ntree = 1000, nodesize = 1, mtry = n^(4/5), nsplit = 10, importance = TRUE)
    obj_error = get.cindex(obj$yvar[,1], obj$yvar[,2], obj$predicted.oob)
    obj_cor = cor(obj$predicted[!(forest_data$time_to_event < 265 & forest_data$event_stat == 0)], forest_data$time_to_event[!(forest_data$time_to_event < 265 & forest_data$event_stat == 0)])
    errors = rbind(errors, c("coxph", n, trial, obj_error, length(list_coxph[[n]][["topvars"]]), obj_cor))
  }
  
  for(trial in seq(1, 5)){
    print(paste("DESeq2", trial))
    obj = rfsrc(Surv(time_to_event, event_stat)~., forest_data[, c(list_deseq2[[n]][["topvars"]], "time_to_event", "event_stat")], ntree = 1000, nodesize = 1, mtry = n^(4/5), nsplit = 10, importance = TRUE)
    obj_error = get.cindex(obj$yvar[,1], obj$yvar[,2], obj$predicted.oob)
    obj_cor = cor(obj$predicted[!(forest_data$time_to_event < 265 & forest_data$event_stat == 0)], forest_data$time_to_event[!(forest_data$time_to_event < 265 & forest_data$event_stat == 0)])
    errors = rbind(errors, c("deseq2", n, trial, obj_error, length(list_deseq2[[n]][["topvars"]]), obj_cor))
  }
  
  for(trial in seq(1, 5)){
    print(paste("Uniform", trial))
    obj = rfsrc(Surv(time_to_event, event_stat)~., forest_data[, c(list_uniform[[n]][["topvars"]], "time_to_event", "event_stat")], ntree = 1000, nodesize = 1, mtry = n^(4/5), nsplit = 10, importance = TRUE)
    obj_error = get.cindex(obj$yvar[,1], obj$yvar[,2], obj$predicted.oob)
    obj_cor = cor(obj$predicted[!(forest_data$time_to_event < 265 & forest_data$event_stat == 0)], forest_data$time_to_event[!(forest_data$time_to_event < 265 & forest_data$event_stat == 0)])
    errors = rbind(errors, c("uniform", n, trial, obj_error, length(list_uniform[[n]][["topvars"]]), obj_cor))
  }
  
  for(trial in seq(1, 5)){
    print(paste("VIMP", trial))
    obj = rfsrc(Surv(time_to_event, event_stat)~., forest_data[, c(list_vimp[[n]][["topvars"]], "time_to_event", "event_stat")], ntree = 1000, nodesize = 1, mtry = n^(4/5), nsplit = 10, importance = TRUE)
    obj_error = get.cindex(obj$yvar[,1], obj$yvar[,2], obj$predicted.oob)
    obj_cor = cor(obj$predicted[!(forest_data$time_to_event < 265 & forest_data$event_stat == 0)], forest_data$time_to_event[!(forest_data$time_to_event < 265 & forest_data$event_stat == 0)])
    errors = rbind(errors, c("vimp", n, trial, obj_error, length(list_vimp[[n]][["topvars"]]), obj_cor))
  }
  
  for(trial in seq(1,5)){
    print(paste("Union", trial))
    obj = rfsrc(Surv(time_to_event, event_stat)~., forest_data[, c(union_list[[n]], "time_to_event", "event_stat")], ntree = 1000, nodesize = 1, mtry = n^(4/5), nsplit = 10, importance = TRUE)
    obj_error = get.cindex(obj$yvar[,1], obj$yvar[,2], obj$predicted.oob)
    obj_cor = cor(obj$predicted[!(forest_data$time_to_event < 265 & forest_data$event_stat == 0)], forest_data$time_to_event[!(forest_data$time_to_event < 265 & forest_data$event_stat == 0)])
    errors = rbind(errors, c("union", n, trial, obj_error, length(union_list[[n]]), obj_cor))
  }
  
  for(trial in seq(1,5)){
    print(paste("Intersection", trial))
    obj = rfsrc(Surv(time_to_event, event_stat)~., forest_data[, c(intersection_list[[n]], "time_to_event", "event_stat")], ntree = 1000, nodesize = 1, mtry = n^(4/5), nsplit = 10, importance = TRUE)
    obj_error = get.cindex(obj$yvar[,1], obj$yvar[,2], obj$predicted.oob)
    obj_cor = cor(obj$predicted[!(forest_data$time_to_event < 265 & forest_data$event_stat == 0)], forest_data$time_to_event[!(forest_data$time_to_event < 265 & forest_data$event_stat == 0)])
    errors = rbind(errors, c("intersection", n, trial, obj_error, length(intersection_list[[n]]), obj_cor))
  }
  rm(obj, obj_cor, obj_error)
}
colnames(errors) = c("method", "genes", "trial", "error", "length", "cor")
errors$genes = as.numeric(errors$genes)
errors$error = as.numeric(errors$error)
errors$length = as.numeric(errors$length)
errors$cor = as.numeric(errors$cor)
errors$method = factor(errors$method, levels = c("coxph", "deseq2", "uniform", "vimp", "intersection", "union"))

#Do some plotting
errors_error <- data_summary(errors, varname="error", groupnames=c("method", "genes")) #Could I come up with better variable names? Screw you
ggplot(errors_error[!errors_error$method %in% c("union", "intersection"), ], aes(x=as.numeric(genes), y=error, group=method, color=method)) + 
  geom_line() +
  geom_point(size = 3/fig_factor)+
  geom_errorbar(aes(ymin=error-sd, ymax=error+sd), width=.3/fig_factor,
                position=position_dodge(0.05)) + theme_linedraw() + theme(axis.title = element_text(size = 16/fig_factor), axis.text = element_text(size = 12/fig_factor), legend.title = element_text(size = 12/fig_factor), legend.text = element_text(size = 10/fig_factor)) + guides(color = guide_legend(title = "Weighting Method")) + scale_color_manual(values = c("indianred1", "royalblue", "forestgreen", "gold"), labels = c("Cox Proportional Hazard", "DESeq2", "Uniform/Random", "VIMP")) + ylab("Forest Error") + xlab("Number of Input Genes")
ggplot(errors_error, aes(x=as.numeric(genes), y=error, group=method, color=method)) + 
  geom_line() +
  geom_point(size = 3/fig_factor)+
  geom_errorbar(aes(ymin=error-sd, ymax=error+sd), width=.3/fig_factor,
                position=position_dodge(0.05)) + theme_linedraw() + theme(axis.title = element_text(size = 16/fig_factor), axis.text = element_text(size = 12/fig_factor), legend.title = element_text(size = 12/fig_factor), legend.text = element_text(size = 10/fig_factor)) + guides(color = guide_legend(title = "Weighting Method")) + scale_color_manual(values = c("indianred1", "royalblue", "forestgreen", "gold", "black", "gray"), labels = c("Cox Proportional Hazard", "DESeq2", "Uniform/Random", "VIMP", "Intersection", "Union")) + ylab("Forest Error") + xlab("Number of Input Genes")

errors_cor <- data_summary(errors, varname="cor", groupnames=c("method", "genes"))
ggplot(errors_cor[!errors_cor$method %in% c("union", "intersection"), ], aes(x=as.numeric(genes), y=-1*cor, group=method, color=method)) + 
  geom_line() +
  geom_point(size = 3/fig_factor)+
  geom_errorbar(aes(ymin=(-1*cor)-sd, ymax=(-1*cor)+sd), width=.3/fig_factor,
                position=position_dodge(0.05)) + theme_linedraw() + theme(axis.title = element_text(size = 16/fig_factor), axis.text = element_text(size = 12/fig_factor), legend.title = element_text(size = 12/fig_factor), legend.text = element_text(size = 10/fig_factor)) + guides(color = guide_legend(title = "Weighting Method")) + scale_color_manual(values = c("indianred1", "royalblue", "forestgreen", "gold"), labels = c("Cox Proportional Hazard", "DESeq2", "Uniform/Random", "VIMP")) + ylab("Correlation with Survival") + xlab("Number of Input Genes")
ggplot(errors_cor, aes(x=as.numeric(genes), y=-1*cor, group=method, color=method)) + 
  geom_line() +
  geom_point(size = 3/fig_factor)+
  geom_errorbar(aes(ymin=(-1*cor)-sd, ymax=(-1*cor)+sd), width=.3/fig_factor,
                position=position_dodge(0.05)) + theme_linedraw() + theme(axis.title = element_text(size = 16/fig_factor), axis.text = element_text(size = 12/fig_factor), legend.title = element_text(size = 12/fig_factor), legend.text = element_text(size = 10/fig_factor)) + guides(color = guide_legend(title = "Weighting Method")) + scale_color_manual(values = c("indianred1", "royalblue", "forestgreen", "gold", "black", "gray"), labels = c("Cox Proportional Hazard", "DESeq2", "Uniform/Random", "VIMP", "Intersection", "Union")) + ylab("Correlation with Survival") + xlab("Number of Input Genes")

#My take-away from the above analysis is that firstly, there is clearly a sweet spot in terms of number of input genes. You need enough to be able to include the most informative genes as part of the variable hunting search space. However, go too high, and the algorithm probably finds too many false positives because it is overwhelmed by noise. I think the optimal number comes anywhere around 3000-4000, so I chose 4000. Secondly, using the intersectino gene list seems to perform better than any individual method in terms of minimizing forest error, while being very comparable in terms of correlation with survival. That speaks to using the intersection list as our final optimized gene list for further analysis.

#Now, above I used the default medium conservativeness settings to generate the gene lists. I was a little curious to see how the liberalized lists did, so I wrote some quick code to check that (just at n = 4000 because that's the only n at which I re-ran the liberal code above). I think this is mostly for entertainment so I commented it out.
#combined = c(liberal[["coxph"]][["topvars"]], liberal[["deseq2"]][["topvars"]], liberal[["uniform"]][["topvars"]], liberal[["vimp"]][["topvars"]])
#tabled = table(combined)
#tabled = table(combined)
#union_liberal_list = unique(combined)
#intersection_liberal_list = names(tabled)[tabled >= 3]
#rm(combined, tabled)

#errors$method = factor(errors$method, levels = c("coxph", "deseq2", "uniform", "vimp", "intersection", "union", "intersection liberal", "union liberal"))
#for(trial in seq(1,5)){
#  print(paste("Intersection Liberal", trial))
#  obj = rfsrc(Surv(time_to_event, event_stat)~., forest_data[, c(intersection_liberal_list, "time_to_event", "event_stat")], ntree = 1000, nodesize = 1, mtry = n^(4/5), nsplit = 10, importance = TRUE)
#  obj_error = get.cindex(obj$yvar[,1], obj$yvar[,2], obj$predicted.oob)
#  obj_cor = cor(obj$predicted[!(forest_data$time_to_event < 265 & forest_data$event_stat == 0)], forest_data$time_to_event[!(forest_data$time_to_event < 265 & forest_data$event_stat == 0)])
# errors = rbind(errors, c("intersection liberal", 4000, trial, obj_error, length(intersection_liberal_list), obj_cor))
#}
#for(trial in seq(1,5)){
# print(paste("Union Liberal", trial))
# obj = rfsrc(Surv(time_to_event, event_stat)~., forest_data[, c(union_liberal_list, "time_to_event", "event_stat")], ntree = 1000, nodesize = 1, mtry = n^(4/5), nsplit = 10, importance = TRUE)
# obj_error = get.cindex(obj$yvar[,1], obj$yvar[,2], obj$predicted.oob)
# obj_cor = cor(obj$predicted[!(forest_data$time_to_event < 265 & forest_data$event_stat == 0)], forest_data$time_to_event[!(forest_data$time_to_event < 265 & forest_data$event_stat == 0)])
# errors = rbind(errors, c("union liberal", 4000, trial, obj_error, length(union_liberal_list), obj_cor))
#}
#rm(obj, obj_cor, obj_error)

#colnames(errors) = c("method", "genes", "trial", "error", "length", "cor")
#errors$genes = as.numeric(errors$genes)
#errors$error = as.numeric(errors$error)
#errors$length = as.numeric(errors$length)
#errors$cor = as.numeric(errors$cor)
#errors$method = factor(errors$method, levels = c("coxph", "deseq2", "uniform", "vimp", "intersection", "union", "intersection liberal", "union liberal"))

#errors_error <- data_summary(errors, varname="error", groupnames=c("method", "genes"))
#ggplot(errors_error, aes(x=as.numeric(genes), y=error, group=method, color=method)) + 
#  geom_line() +
#  geom_point(size = 3/fig_factor)+
#  geom_errorbar(aes(ymin=error-sd, ymax=error+sd), width=.3/fig_factor,
#                position=position_dodge(0.05)) + theme_linedraw() + theme(axis.title = element_text(size = 16/fig_factor), axis.text = element_text(size = 12/fig_factor), legend.title = element_text(size = 12/fig_factor), legend.text = element_text(size = 10/fig_factor)) + guides(color = guide_legend(title = "Weighting Method")) + scale_color_manual(values = c("indianred1", "royalblue", "forestgreen", "gold", "black", "gray", "pink", "orange"), labels = c("Cox Proportional Hazard", "DESeq2", "Uniform/Random", "VIMP", "Intersection", "Union", "Intersection Liberal", "Union Liberal")) + ylab("Forest Error") + xlab("Number of Input Genes")

#errors_cor <- data_summary(errors, varname="cor", groupnames=c("method", "genes"))
#ggplot(errors_cor, aes(x=as.numeric(genes), y=-1*cor, group=method, color=method)) + 
#  geom_line() +
#  geom_point(size = 3/fig_factor)+
#  geom_errorbar(aes(ymin=(-1*cor)-sd, ymax=(-1*cor)+sd), width=.3/fig_factor,
#                position=position_dodge(0.05)) + theme_linedraw() + theme(axis.title = element_text(size = 16/fig_factor), axis.text = element_text(size = 12/fig_factor), legend.title = element_text(size = 12/fig_factor), legend.text = element_text(size = 10/fig_factor)) + guides(color = guide_legend(title = "Weighting Method")) + scale_color_manual(values = c("indianred1", "royalblue", "forestgreen", "gold", "black", "gray", "pink", "orange"), labels = c("Cox Proportional Hazard", "DESeq2", "Uniform/Random", "VIMP", "Intersection", "Union", "Intersectional Liberal", "Union Liberal")) + ylab("Correlation with Survival") + xlab("Number of Input Genes")

###Selecting the optimized gene list###
#Based on the logic above, we can now select our final optimized gene list. This will be the list created as an intersection of the variable lists created from variable hunting with 4000 input genes (sorted by normalized gene variance), across the four weighting methods, with the criterion that the gene had to have appeared in the list generated by at least three of the methods. Because I am not creative, I am labeling this list "best" - yes, I'm aware that I should probably do better.
best = intersection_list[[4000]]

#Creating a PCA plot based on the best genes
pca_hfpef = prcomp(clean_data[best, ])
pca_hfpef_plot = as.data.frame(pca_hfpef$rotation)
pca_hfpef_plot$category = pheno_hfpef$category
ggplot(pca_hfpef_plot, aes(x = PC1, y = PC2, color = category)) + geom_point(size = 4/fig_factor) + theme_linedraw() + theme(axis.title = element_text(size = 16/fig_factor), axis.text = element_text(size = 12/fig_factor), legend.title = element_text(size = 12/fig_factor), legend.text = element_text(size = 10/fig_factor)) + guides(color = guide_legend(title = "Year of Event")) + scale_color_manual(values = c("indianred1", "royalblue", "orange", "gold")) + ylab("PC2") + xlab("PC1")

#Creating a plot of predicted survival score (based on the best genes) vs the actual time to event
obj = rfsrc(Surv(time_to_event, event_stat)~., forest_data[, c(best, "time_to_event", "event_stat")], ntree = 1000, nodesize = 1, mtry = n^(4/5), nsplit = 10, importance = TRUE)
pheno_hfpef$predicted = obj$predicted
gene_cor = numeric()
for(gene in best){
  gene_cor[gene] = cor(clean_data[gene, ], obj$predicted)
}
ggplot(pheno_hfpef, aes(x = time_to_event, y = predicted, color = category)) + geom_point(size = 4/fig_factor) + theme_linedraw() + theme(axis.title = element_text(size = 16/fig_factor), axis.text = element_text(size = 12/fig_factor), legend.title = element_text(size = 12/fig_factor), legend.text = element_text(size = 10/fig_factor)) + guides(color = guide_legend(title = "Year of Event")) + scale_color_manual(values = c("indianred1", "royalblue", "orange", "gold")) + ylab("Predicted Survival Score") + xlab("Time to Event")

#Heatmap of genes with the patients sorted by their time to event
hfpef_heatmap = pheatmap(clean_data[best, order(pheno_hfpef$predicted, decreasing = TRUE)], cluster_rows = TRUE, cluster_cols = FALSE, show_rownames = TRUE, show_colnames = FALSE, annotation_col = pheno_hfpef[order(pheno_hfpef$predicted, decreasing = TRUE), c("category"), drop = FALSE], scale = "row", color = rev(colorRampPalette(brewer.pal(11, "RdBu"))(100)), fontsize = 10/fig_factor)

#We can plot these genes for the control, HFpEF, and HFrEF samples just for fun
rld_test = assay(rld)
rownames(rld_test) = rownames(clean_data)
pheno$disease = factor(pheno$disease, levels = c("HFrEF", "HFpEF", "Normal"))
pheatmap(rld_test[intersection_list[[4000]][hfpef_heatmap$tree_row$order], rownames(pheno)[pheno$disease %in% c("HFrEF", "Normal") & pheno$tissue == "RVS"]], cluster_rows = FALSE, cluster_cols = FALSE, show_rownames = TRUE, show_colnames = FALSE, annotation_col = pheno[rownames(pheno)[pheno$disease %in% c("HFrEF", "Normal") & pheno$tissue == "RVS"], "disease", drop = FALSE], scale = "row", color = rev(colorRampPalette(brewer.pal(11, "RdBu"))(100)))
fold_changes = data.frame(HFpEF_survival = survival[rownames(clean_data), ]$log2FoldChange, HFrEF = hfref[rownames(clean_data), ]$log2FoldChange)
rownames(fold_changes) = rownames(clean_data)

###Loading in clinical parameters###
clinical_data = read.csv("~/Documents/Research/HFpEF/clin_param.csv", as.is = TRUE)
rownames(clinical_data) = paste("s", clinical_data$SequenceRunId, sep = "")
clinical_data = clinical_data[rownames(pheno_hfpef), ]
clinical_data = clinical_data[, 7:115]
clinical_data$event_stat = as.numeric(pheno_hfpef$event_stat)
clinical_data$time_to_event = pheno_hfpef$time_to_event
rownames_save = rownames(clinical_data)
clinical_data = as.data.frame(unclass(clinical_data),stringsAsFactors=TRUE)
rownames(clinical_data) = rownames_save
rm(rownames_save)

#We can now run correlations between all of the clinical parameters and the genes - of course, not all of these will be totally meaningful, but it's quick to calculate so we can run this now and use as needed later
var_genes = names(sort(testVars, decreasing = TRUE))[1:4000]
clin_param_marker_cor = expand.grid(var_genes, colnames(clinical_data))
clin_param_marker_cor$cor = 0
for(i in seq(1:nrow(clin_param_marker_cor))){
  clin_param_marker_cor[i, ]$cor = cor(clean_data[clin_param_marker_cor[i, 1], ], as.numeric(clinical_data[, clin_param_marker_cor[i, 2]]), method = "spearman", use="complete.obs")
}  

#For further analysis, I narrowed down the clinical parameters to a short list focused on hemodynamic and echo parameters, though with some other values in there. I also tried to remove parameters that were either duplicates or likely to overlap (e.g. Cr and eGFR), though in some cases I kept in duplicates (for example, there are several pro-BNP measurements, as well as NPPB).
best_clinical_parameters = read.table("best_clin_param.txt", as.is = TRUE)$V1
best_clinical_data = cbind(clinical_data[, best_clinical_parameters], t(clean_data[best, ]))
best_clinical_data$time_to_event = forest_data$time_to_event
best_clinical_data$event_stat = forest_data$event_stat

#Now, we can create random forests from combinations of the optimal gene list and/or clinical parameters as above to compare their performance.
#Forest with just the genes
obj = rfsrc(Surv(time_to_event, event_stat)~., best_clinical_data[, colnames(best_clinical_data) %in% c(best, "event_stat", "time_to_event")], importance = TRUE, ntree = 1000, nodesize = 1, mtry = length(best)^(4/5), nsplit = 10)

#Forest with a random sampling of 33 genes not in the optimal gene list (as a control)
obj = rfsrc(Surv(time_to_event, event_stat)~., forest_data[, c(sample(var_genes[!var_genes %in% best], size = 33), "event_stat", "time_to_event")], importance = TRUE, ntree = 1000, nodesize = 1, mtry = 33^(4/5), nsplit = 10)

#Forest with just the clinical parameters
obj = rfsrc(Surv(time_to_event, event_stat)~., best_clinical_data[, !colnames(best_clinical_data) %in% best], importance = TRUE, ntree = 1000, nodesize = 1, mtry = length(best_clinical_parameters)^(4/5), nsplit = 10)

#Forest combining the genes and clinical parameters
obj = rfsrc(Surv(time_to_event, event_stat)~., best_clinical_data, importance = TRUE, ntree = 1000, nodesize = 1, mtry = (length(best) + length(best_clinical_parameters))^(4/5), nsplit = 10)

#Forest combining the clinical parameters with a set of 33 randomly selected genes NOT in our optimal gene list (as a sort of control group)
obj = rfsrc(Surv(time_to_event, event_stat)~., cbind(best_clinical_data[, !colnames(best_clinical_data) %in% best], forest_data[, sample(var_genes[!var_genes %in% best], size = 33)]), importance = TRUE, ntree = 1000, nodesize = 1, mtry = (length(best) + length(best_clinical_parameters))^(4/5), nsplit = 10)

#We can do each of the above for multiple trials and combine to make one nice figure. 
param_cindex = data.frame(matrix(ncol = 3))
for(trial in seq(1, 10)){
  obj = rfsrc(Surv(time_to_event, event_stat)~., best_clinical_data[, colnames(best_clinical_data) %in% c(best, "event_stat", "time_to_event")], importance = TRUE, ntree = 1000, nodesize = 1, mtry = length(best)^(4/5), nsplit = 10)
  param_cindex = rbind(param_cindex, c("best_genes", trial, 1-get.cindex(obj$yvar[,1], obj$yvar[,2], obj$predicted.oob)))
  print(paste("Best Genes", trial))
}
for(trial in seq(1, 10)){
  obj = rfsrc(Surv(time_to_event, event_stat)~., forest_data[, c(sample(var_genes[!var_genes %in% best], size = 33), "event_stat", "time_to_event")], importance = TRUE, ntree = 1000, nodesize = 1, mtry = 33^(4/5), nsplit = 10)
  param_cindex = rbind(param_cindex, c("random_genes", trial, 1-get.cindex(obj$yvar[,1], obj$yvar[,2], obj$predicted.oob)))
  print(paste("Random Genes", trial))
}
for(trial in seq(1, 10)){
  obj = rfsrc(Surv(time_to_event, event_stat)~., best_clinical_data[, !colnames(best_clinical_data) %in% best], importance = TRUE, ntree = 1000, nodesize = 1, mtry = length(best_clinical_parameters)^(4/5), nsplit = 10)
  param_cindex = rbind(param_cindex, c("clin_param", trial, 1-get.cindex(obj$yvar[,1], obj$yvar[,2], obj$predicted.oob)))
  print(paste("Clinical Parameters", trial))
}
for(trial in seq(1, 10)){
  obj = rfsrc(Surv(time_to_event, event_stat)~., best_clinical_data, importance = TRUE, ntree = 1000, nodesize = 1, mtry = (length(best) + length(best_clinical_parameters))^(4/5), nsplit = 10)
  param_cindex = rbind(param_cindex, c("clin_param_best_genes", trial, 1-get.cindex(obj$yvar[,1], obj$yvar[,2], obj$predicted.oob)))
  print(paste("Clinical Parameters + Best Genes", trial))
}
for(trial in seq(1, 10)){
  obj = rfsrc(Surv(time_to_event, event_stat)~., cbind(best_clinical_data[, !colnames(best_clinical_data) %in% best], forest_data[, sample(var_genes[!var_genes %in% best], size = 33)]), importance = TRUE, ntree = 1000, nodesize = 1, mtry = (length(best) + length(best_clinical_parameters))^(4/5), nsplit = 10)
  param_cindex = rbind(param_cindex, c("clin_param_random_genes", trial, 1-get.cindex(obj$yvar[,1], obj$yvar[,2], obj$predicted.oob)))
  print(paste("Clinical Parameters + Random Genes", trial))
}
colnames(param_cindex) = c("parameter", "trial", "c_index")
param_cindex = param_cindex[2:51, ]
param_cindex$c_index = as.numeric(param_cindex$c_index)
param_cindex$parameter = factor(param_cindex$parameter, levels = c("best_genes", "random_genes", "clin_param", "clin_param_best_genes", "clin_param_random_genes"))
param_cindex$pretty_names = c(rep("Optimal Genes", 10), rep("Random Genes", 10), rep("Clinical Parameters", 10), rep("Optimal Genes + Clinical Parameters", 10), rep("Random Genes + Clinical Parameters", 10))
ggplot(param_cindex, aes(x = pretty_names, y = c_index)) + geom_boxplot(linewidth = 1/fig_factor) + ylim(0, 1.0) + theme_linedraw() + theme(axis.title = element_text(size = 16/fig_factor), axis.text = element_text(size = 12/fig_factor), axis.text.x = element_text(angle = 45), legend.title = element_text(size = 12/fig_factor), legend.text = element_text(size = 10/fig_factor)) + ylab("C-index") + xlab("Included Parameters") + scale_x_discrete(labels = function(x) str_wrap(x, width = 20))

#As another neat test, we can create forests with the clinical parameters but spiking in the genes from the optimal gene list - this will pretty clearly show that for each gene we add, our predictive capacity improves pretty dramatically
obj = rfsrc(Surv(time_to_event, event_stat)~., best_clinical_data[, colnames(best_clinical_data) %in% c(best, "event_stat", "time_to_event")], importance = TRUE, ntree = 1000, nodesize = 1, mtry = length(best)^(4/5), nsplit = 10)
imp = sort(obj$importance, decreasing = TRUE) #We start by ordering our genes by importance from a test forest
param_cindex_add = data.frame(matrix(ncol = 3))
for(gene_no in seq(0, length(best))){
  print(gene_no)
  for(trial in seq(1, 10)){
    obj = rfsrc(Surv(time_to_event, event_stat)~., best_clinical_data[c(best_clinical_parameters, names(imp)[0:gene_no], "time_to_event", "event_stat")], importance = TRUE, ntree = 1000, nodesize = 1, mtry = (gene_no + length(best_clinical_parameters))^(4/5), nsplit = 10)
    param_cindex_add = rbind(param_cindex_add, c(gene_no, trial, 1-get.cindex(obj$yvar[,1], obj$yvar[,2], obj$predicted.oob)))
  }
}
param_cindex_add = param_cindex_add[2:nrow(param_cindex_add), ]
colnames(param_cindex_add) = c("genes", "trial", "c_index")
param_cindex_add$c_index = as.numeric(param_cindex_add$c_index)
ggplot(param_cindex_add, aes(x = as.factor(genes), y = c_index)) + geom_boxplot(linewidth = 0.5/fig_factor) + ylim(0, 1.0) + theme_linedraw() + theme(axis.title = element_text(size = 16/fig_factor), axis.text = element_text(size = 12/fig_factor), legend.title = element_text(size = 12/fig_factor), legend.text = element_text(size = 10/fig_factor)) + ylab("C-index") + xlab("Added Number of Genes") + scale_x_discrete(labels = function(x) str_wrap(x, width = 20))

#We can also compare the importance of individual variables - particularly, comparing the clinical variables to the identified genes. Here, I created a forest using both the best genes and the clinical parameters, 
obj = rfsrc(Surv(time_to_event, event_stat)~., best_clinical_data, importance = TRUE, ntree = 1000, nodesize = 1, mtry = (length(best) + length(best_clinical_parameters))^(4/5), nsplit = 10)
obj2 = max.subtree(obj)
importance = data.frame(row.names = obj$xvar.names, importance = obj$importance, minimal_depth = obj2$order[, 1])
importance$variable = "Clinical Parameter"
importance[rownames(importance) %in% best, ]$variable = "Gene"
importance$names = ""
importance[rownames(importance) %in% c("ADAMTSL2", "DNAJB9", "KATNAL1", "ADRB1", "METRNL", "BMP6", "LVEDD..cm", "PAmean..mmHg", "CI..L.min.m2", "NTproBNP..pg.mL", "NYHA.Class", "Age", "eGFR..mL.min.1.73m2", "Hemoglobin"), ]$names = c("Age", "NYHA.Class", "LVEDD", "eGFR", "NT-proBNP", "PA Mean", "CI", "Hgb", "ADAMTSL2", "ADRB1", "BMP6", "DNAJB9", "KATNAL1", "METRNL")
importance = importance[order(importance$importance, decreasing = TRUE), ]
ggplot(importance, aes(x = seq(1, length(best) + length(best_clinical_parameters)), y = importance, color = variable)) + geom_point() + geom_label_repel(aes(label = names), size = 3, box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50', max.overlaps = 10, nudge_x = 2, nudge_y = 0.015) + theme_bw() + ylab("VIMP") + xlab("Variable")
importance = importance[order(importance$minimal_depth), ]
ggplot(importance, aes(x = seq(1, length(best) + length(best_clinical_parameters)), y = minimal_depth, color = variable)) + geom_point() + geom_label_repel(aes(label = names), size = 3, box.padding = 0.35, point.padding = 0.5, segment.color = 'grey50', max.overlaps = 10, nudge_x = 2, nudge_y = 0.5) + theme_bw() + ylab("Minimal Depth") + xlab("Variable")

###Survival analysis using optimized gene list###
#Now what we will do is to use the identified optimized gene list (33 genes) to split the patients into groups. We can do this in an unbiased manner as follows. Each gene is either positively or negatively correlated to some degree with survival. Suppose we scale each gene from 0 to 1 to correlate to survival. We can then classify patients based on a combined gene score. Based on this approach, I split them into three groups, and then did survival analysis.

best_gene_scale = clean_data[best, ]
best_gene_scale[gene_cor > 0, ] = best_gene_scale[gene_cor > 0, ]*-1
best_gene_rank = rowRanks(best_gene_scale)
gene_score = colMeans(best_gene_rank)
pheno_hfpef$gene_score = gene_score
gene_score_interval = (max(gene_score) - min(gene_score))/3
pheno_hfpef$gene_score_group = 2
pheno_hfpef[pheno_hfpef$gene_score < (min(gene_score) + gene_score_interval), ]$gene_score_group = 1
pheno_hfpef[pheno_hfpef$gene_score > (max(gene_score) - gene_score_interval), ]$gene_score_group = 3
#pheno_hfpef$gene_score_group = factor(pheno_hfpef$gene_score_group, levels = c(3, 2, 1))
pheno_hfpef$gene_score_group = factor(pheno_hfpef$gene_score_group, levels = c(1, 2, 3))
fit <- survfit(Surv(time_to_event, event_stat) ~ gene_score_group,
               data = pheno_hfpef)
ggsurvplot(fit, data = pheno_hfpef, risk.table = FALSE, ggtheme = theme_bw(), conf.int = TRUE)
pheno_hfpef$renal_function = clinical_data$eGFR..mL.min.1.73m2
coxph(Surv(time_to_event, event_stat) ~ gene_score_group + sex + age + renal_function,
      data = pheno_hfpef) #This allows us to test for significance using Cox while account for sex, age, and renal function (in the form of eGFR) as confounders

###Single cell RNA-seq Analysis###
#We now do some exploratory analysis of our identified survival genes in myocardial single cell RNA-seq. At the time I was putting this code together, Dr. Hahn's single cell data of HFpEF was not publicly available, and so this analysis is based on Litvinukova et al. (2020). You can readily adapt this code to be used for other datasets too, however; and as of the finalized version of this code, Dr. Hahn's data is available through the preprint on bioRxiv.

#Note - the files provided by the authors are in h5ad format, primarily for use through Python. Since I'm working through R here (too lazy to go back and forth between languages), my first step was to convert the file to Seurat output using sceasy. 
sceasy::convertFormat("~/Downloads/Global_raw.h5ad", from="anndata", to="seurat", outFile='~/Downloads/Global_raw.rds')
data = readRDS("~/Downloads/Global_raw.rds")

#Next, some clean-up of the metadata. This was largely done to synchronize this dataset with Dr. Hahn's, and to make it a little easier to read
data@meta.data$cell_type = factor(data@meta.data$cell_type, levels= levels(data@meta.data$cell_type)[c(9, 1, 2, 4, 3, 11, 5, 12, 7, 8, 6, 10)])
data@meta.data$good_celltype = as.character(data@meta.data$cell_type)
data@meta.data[data@meta.data$good_celltype %in% c("Atrial Cardiomyocyte", "Ventricular Cardiomyocyte"), ]$good_celltype = "Cardiomyocyte"
data@meta.data[data@meta.data$good_celltype %in% c("Lymphatic Endothelial cell"), ]$good_celltype = "Lymphatic cell"
data@meta.data$good_celltype = factor(data@meta.data$good_celltype, levels = c("Adipocyte", "Cardiomyocyte", "Mesothelial cell", "Fibroblast", "Endothelial cell", "Lymphoid", "Lymphatic cell", "Mast cell", "Neural cell", "Mural cell", "Myeloid"))

#The Litvinukova dataset uses ensembl gene names, so I converted quickly
gene_ens = c("ENSG00000102781", "ENSG00000126790", "ENSG00000173599", "ENSG00000277363", "ENSG00000106772", "ENSG00000136021", "ENSG00000167969", "ENSG00000113389", "ENSG00000138029", "ENSG00000168952", "ENSG00000006468", "ENSG00000006062", "ENSG00000043591", "ENSG00000146926", "ENSG00000141639", "ENSG00000134986", "ENSG00000236751", "ENSG00000185028", "ENSG00000040199", "ENSG00000128590", "ENSG00000172201", "ENSG00000162738", "ENSG00000197859", "ENSG00000164877", "ENSG00000000971", "ENSG00000167434", "ENSG00000153162", "ENSG00000014914", "ENSG00000110092", "ENSG00000164318", "ENSG00000244405", "ENSG00000176845", "ENSG00000198585")
genes2 =  c("KATNAL1", "L3HYPDH", "PC", "SRCIN1", "PRUNE2", "SCYL2", "ECI1", "NPR3", "HADHB", "STXBP6", "ETV1", "MAP3K14", "ADRB1", "ASB10", "MAPK4", "NREP",  "LINC01186", "LRRC14B",  "PHLPP2", "DNAJB9", "ID4", "VANGL2", "ADAMTSL2", "MICALL2", "CFH","CA4", "BMP6", "MTMR11",   "CCND1",   "EGFLAM", "ETV5", "METRNL", "NUDT16")
names(gene_ens) = genes2
DotPlot(data, group.by = "good_celltype", features = gene_ens, dot.scale = 6/fig_factor, scale.max = 30, col.max = 5) + RotatedAxis() + theme(axis.text.x = element_text(size = 12/fig_factor), axis.text.y = element_text(size = 16/fig_factor), legend.position = "none") + scale_x_discrete(labels = genes2)

###Cross-validation with other datasets###
#Cross-validation is extremely challenging because there just aren't that many datasets of human HFpEF tissue, and certainly none that have paired clinical data. As a sort of pseudo-cross validation, I used the data from Das et al. (2019). Keep in mind - this was a very inexact comparison. That dataset used patients who were going for elective CABG, and retrospectively labeled 5 of them HFpEF-proxy based on echo and proBNP. There was no data to indicate whether these patients were decompensated. My rough assumption was that the HFpEF-proxy patients should have worse predicted survival scores than the controls. Again - this is very inexact but the best available cross-comparison.

#Cleans the Das data and our data; we log2scale their data to make it comparable to ours, and then scale both datasets so that the numbers are a bit more comparable to each other (otherwise the absolute counts would not make sense)
scale_data = as.data.frame(scale(forest_data[, best]))
das_data = read.table("das_data.txt", as.is = TRUE, header = TRUE, row.names = 1)
das_scale_data = as.data.frame(scale(t(log2(das_data[gene_ens[best], ]))))
colnames(das_scale_data) = best
das_scale_data = das_scale_data[, !is.na(das_scale_data[1, ])]
scale_data = scale_data[, colnames(scale_data) %in% colnames(das_scale_data)]
scale_data$time_to_event = forest_data$time_to_event
scale_data$event_stat = forest_data$event_stat
obj = rfsrc(Surv(time_to_event, event_stat)~., scale_data, importance = TRUE, ntree = 1000, nodesize = 1, mtry = length(best)^(4/5), nsplit = 10)
das_prediction = predict(obj, das_scale_data)

surv_func = as.data.frame(das_prediction$survival)
colnames(surv_func) = das_prediction$time.interest
surv_func$`1095` = surv_func$`619` 
surv_time = numeric(length = nrow(surv_func))
for(i in seq(1, nrow(surv_func))){
  surv_time[i] = integrate(approxfun(as.numeric(colnames(surv_func)), surv_func[i, ]), lower = 0, upper = max(as.numeric(colnames(surv_func))))$value
} #This basically integrates the area under the survival curve to put the predicted survival on an order of days, comparable to our actual input data
das_pheno = data.frame(row.names = colnames(das_data), group = sapply(str_split(colnames(das_data), "_"), "[[", 3), survival = surv_time)
ggplot(das_pheno, aes(x = group, y = survival)) + geom_boxplot(linewidth = 1/fig_factor) + geom_jitter(size = 5/fig_factor) + theme_linedraw() + theme(axis.title = element_text(size = 16/fig_factor), axis.text = element_text(size = 12/fig_factor), legend.title = element_text(size = 12/fig_factor), legend.text = element_text(size = 10/fig_factor)) + ylab("Predicted Survival") + xlab("Group") 