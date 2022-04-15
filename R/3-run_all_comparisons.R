# --------------------------------------------------------------------------- #
# Script 3: Comparisons
#
# Compare +/- stretch, EED vs Healthy, etc
# --------------------------------------------------------------------------- #

library(limma)
library(plyr)
library(RColorBrewer)
PLOT_COLS = brewer.pal(11, 'Spectral')

# Volcano plot function
volcano_plot <- function(log2_fc, p_adj, 
                         log2_fc_threshold = 0.5, p_adj_threshold = 0.05, 
                         plot_main = 'Volcano', ymax = 25){
  
  if (length(log2_fc) != length(p_adj)){
    return('Log2FC and pAdj vectors are different lengths!')
  } else {
    
    # Red, teal, grey
    plot_cols = c(PLOT_COLS[2], '#27adde', '#C3C5C7')
    cols = rep(plot_cols[3], length(p_adj))
    cols[log2_fc > log2_fc_threshold & p_adj < p_adj_threshold] = plot_cols[1]
    cols[log2_fc < -log2_fc_threshold & p_adj < p_adj_threshold] = plot_cols[2]
    
    if (is.null(ymax)){
      ymax = max(-log10(p_adj), na.rm = T) + 10
    }
    
    plot(log2_fc, -log10(p_adj),
         main = plot_main,
         xlab = 'log2 fold change', ylab = '-log10 FDR-adjusted p-value',
         ylim = c(0, ymax),
         las = 1, pch = 21,
         bg = cols)
    abline(h = -log10(p_adj_threshold),
           col = plot_cols[3])
    abline(v = log2_fc_threshold, col = plot_cols[1])
    abline(v = -log2_fc_threshold, col = plot_cols[2])
  }
}

# Read sample data
sample_annotations = read.csv('Data/Cleaned/All_sample_data.csv', stringsAsFactors = F)
row.names(sample_annotations) = sample_annotations$Sample_Name

# Sort by chip #
sample_annotations = sample_annotations[order(sample_annotations$Chip, decreasing = F),]

# Read metabolite data
metabolite_data = read.csv('Data/Cleaned/All_metabolite_data_origscale_tic_normalized_rescaled.csv',
                           stringsAsFactors = F, row.names = 1)
metabolite_annotations = read.csv('Data/Cleaned/All_metabolite_annotations.csv', stringsAsFactors = F)

# Differential expression comparisons to perform 
comparisons = data.frame(Comparison_ID = 1:4,
                         Experimental_Group_Name = c('EED_stretch', 'EED_nostretch', 'EED', 'EED Basal'),
                         Control_Group_Name = c('Healthy_stretch', 'Healthy_nostretch',
                                                'Healthy', 'Healthy Basal'),
                         Experimental_Group = c('3a_3b', '4a_4b', '3a_3b_4a_4b', '3b_4b'),
                         Control_Group = c('9a_9b', '10a_10b', '9a_9b_10a_10b', '9b_10b'),
                         stringsAsFactors = F)

sample_data = data.frame(Sample_ID = names(metabolite_data),
                         Group = sample_annotations[names(metabolite_data),'Group_Number'],
                         stringsAsFactors = F)

# PARAMS for analysis
LOG2_FC_THRESHOLD = 0.15 # Approx = 10% difference in abundance
P_ADJ_THRESHOLD = 0.05

for (i in 1:nrow(comparisons)){
  
  CTRL_GROUP = comparisons[i, 'Control_Group']
  EXPT_GROUP = comparisons[i, 'Experimental_Group']
  cat('Comparing', EXPT_GROUP, 'versus', CTRL_GROUP, '\n')
  
  if (grepl('_', CTRL_GROUP)){
    all_samples = c(strsplit(CTRL_GROUP, '_')[[1]],
                    strsplit(EXPT_GROUP, '_')[[1]])
    sample_data_tmp = sample_data[sample_data$Group %in% all_samples,]
    sample_data_tmp$Group = ifelse(sample_data_tmp$Group %in% strsplit(CTRL_GROUP, '_')[[1]], CTRL_GROUP, EXPT_GROUP)
    
  }else{
    sample_data_tmp = sample_data[grepl(CTRL_GROUP, sample_data$Group) | grepl(EXPT_GROUP, sample_data$Group),]
  }
  expression_data_tmp = metabolite_data[,sample_data_tmp$Sample_ID]
  
  design = model.matrix(~0 + as.factor(Group),
                        data = sample_data_tmp)
  
  colnames(design)[grepl(CTRL_GROUP, colnames(design))] = make.names(CTRL_GROUP)
  colnames(design)[grepl(EXPT_GROUP, colnames(design))] = make.names(EXPT_GROUP)
  
  contrast = makeContrasts(get(make.names(EXPT_GROUP))-get(make.names(CTRL_GROUP)), 
                           levels = design)
  colnames(contrast) = paste0(EXPT_GROUP, '_vs_', CTRL_GROUP)
  
  # Make table for diff abundance results
  fit = eBayes(contrasts.fit(lmFit(expression_data_tmp, design), contrast))
  sigtable = topTable(fit, adjust = 'fdr', number = nrow(expression_data_tmp)+1,
                      confint = T)
  sigtable = cbind(Metabolon_Metabolite_ID = row.names(sigtable),
                   sigtable)
  
  # Calculate unmoderated t-test
  fit_unmod = contrasts.fit(lmFit(expression_data_tmp, design), contrast)
  unmod_t = (fit_unmod$coefficients / fit_unmod$stdev.unscaled / fit_unmod$sigma)[,1]
  unmod_df = data.frame(Unmodified_t = unmod_t,
                        Unmodified_pvalue = 2*pt(abs(unmod_t), length(unmod_t)-1, lower=F),
                        Unmodified_DF = fit_unmod$df.residual)
  unmod_df$Metabolon_Metabolite_ID = row.names(unmod_df)
  unmod_df$Unmodified_fdr = p.adjust(unmod_df$Unmodified_pvalue, method = 'fdr')
  sigtable = join(sigtable, unmod_df)
  names(sigtable)[7:8] = c('Modified_t_pvalue', 'Modified_t_fdr')
  
  # Save table
  sigtable = join(metabolite_annotations, sigtable)
  sigtable = sigtable[order(sigtable$Unmodified_fdr),]
  write.csv(sigtable, 
            paste0('Results/Tables/Differential_expression/', EXPT_GROUP, '_vs_', CTRL_GROUP, '.csv'),
            row.names = F, na = '')
  
  # Save volcano plot
  pdf(paste0('Results/Plots/Differential_expression/Volcano_',
             EXPT_GROUP, '_vs_', CTRL_GROUP, '.pdf'),
      height = 7, width = 5)
  volcano_plot(sigtable$logFC, sigtable$Unmodified_fdr, 
               plot_main = paste0('Differentially abundant metabolites\n', EXPT_GROUP, ' vs ', CTRL_GROUP),
               log2_fc_threshold = LOG2_FC_THRESHOLD,
               p_adj_threshold = P_ADJ_THRESHOLD,
               ymax = max(-log10(sigtable$Unmodified_fdr), na.rm = T))
  dev.off()
}