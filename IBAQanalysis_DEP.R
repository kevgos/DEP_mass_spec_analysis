#load packages DEP and dplyr


# choose the proteinGroups.txt file from MaxQuant output
data <- read.delim(file.choose())

#Looking at the dimensions
dim(data)

#Looking at the column names, ensure same as experimental plan for samples
colnames(data)

#Filter out the reverse hits and probable contaminants
data <- data[data$Reverse != "+" & data$Potential.contaminant != "+",]

#Can see some hits were filtered out
dim(data)

#get rid of any duplicated genes

#check if there are duplicated genes
data$Gene.names %>% duplicated() %>% any()

#remove
data %>% group_by(Gene.names) %>% summarize(frequency = n()) %>% 
  arrange(desc(frequency)) %>% filter(frequency > 1)
data_unique <- make_unique(data, "Gene.names", "Protein.IDs", delim = ";")

#check again
data$name %>% duplicated() %>% any()


# generate explerimental plan using table
expPlan <- read.delim(file.choose())
dim(expPlan)
expPlan
typeof(expPlan)


# generate the summarized experiment.

iBAQ_columns <- grep("iBAQ.", colnames(data_unique))

experimental_design <- expPlan

data_se <- make_se(data_unique, iBAQ_columns, experimental_design)

#have a look at it
data_se

# Plot a barplot of the protein identification overlap between samples

# plotting functions also provided

plot_frequency(data_se)

data_filt2 <- filter_missval(data_se, thr = 1)

plot_numbers(data_filt2)

plot_coverage(data_filt2)

data_norm <- normalize_vsn(data_filt2)

plot_normalization(data_filt2, data_norm)

# Plot a heatmap of proteins with missing values
plot_missval(data_filt2)

# Plot intensity distributions and cumulative fraction of proteins with and without missing values
plot_detect(data_filt2)

# Can impute from various given distributions, can get parameters using different methods

impute(data_norm, fun = "")

data_imp <- impute(data_norm, fun = "MinProb", q = 0.01)
plot_imputation(data_norm, data_imp)

# Differential enrichment analysis
# Test manually defined comparisons for samples

data_diff_manual <- test_diff(data_imp, type = "manual", 
                              test = c("sample1", "ctrl"))

# Denote significant proteins based on user defined cutoffs

dep <- add_rejections(data_diff_manual, alpha = 0.5, lfc = log2(0.37))

# Plot the first and second principal components

plot_pca(dep, x = 1, y = 2, n = 500, point_size = 4)

#Can do this for whatever want to plot
plot_volcano(dep, contrast = "sample1_vs_ctrl", label_size = 2, add_names = TRUE)

# generate a results file

data_results <- get_results(dep)

# write to csv file
write.csv(data_results, "results.csv",)

