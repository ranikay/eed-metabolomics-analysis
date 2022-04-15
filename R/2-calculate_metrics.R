# --------------------------------------------------------------------------- #
# Script 2: Calculations
#
# Calculate natural degradation, cellular utilization, and uptake and transfer
# --------------------------------------------------------------------------- #

library(plyr)
library(RColorBrewer)
library(beeswarm)
source('R/heatmap.3.R')
PLOT_COLS = brewer.pal(11, 'Spectral')

# Read metabolite data (origscale data from Metabolon, normalized to total ion count, rescaled)
data_origscale_norm = read.csv('Data/All_metabolite_data_origscale_tic_normalized_rescaled.csv', 
                               stringsAsFactors = F, row.names = 1)

# Read sample data
sample_annotations = read.csv('Data/All_sample_data.csv', stringsAsFactors = F)
row.names(sample_annotations) = sample_annotations$Sample_Name

# Sort by chip #
sample_annotations = sample_annotations[order(sample_annotations$Chip, decreasing = F),]

# Loading control samples
controls = sample_annotations[grepl('Loading media', sample_annotations$Condition), 'Sample_Name']

# Read metabolite annotations
metabolite_annotations = read.csv('Data/All_metabolite_annotations.csv',
                                  stringsAsFactors = F)

# Extract all the different sample groups
for (group_number in unique(sample_annotations$Group_Number)){
  cat(group_number, '\n')
  assign(paste0('group_', group_number),
         sample_annotations[sample_annotations$Group_Number == group_number, 'Sample_Name'])
  
  if(all(get(paste0('group_', group_number)) %in% names(data_origscale_norm))){
    assign(paste0('group_', group_number, '_avg'),
           apply(data_origscale_norm[,get(paste0('group_', group_number))], 1, mean, na.rm = T))
  }
}

na0 <- function(x){
  x[is.na(x)] = 0
  return(x)
}

### Calculation 1: Natural Degradation (D)
# Loading DM-Apical Inlet
# D1=5a-7a
# D2=5a-8a
# Dav=average D1:D2

# There aren't chip #s for the loading media controls so we just use the average here

# D1 = avg(5) - avg(7a)
D = data.frame(Metabolon_Metabolite_ID = row.names(data_origscale_norm),
               D1 = apply(data_origscale_norm[,group_5a], 1, mean, na.rm = T) - 
                 apply(data_origscale_norm[,group_7a], 1, mean, na.rm = T),
               D2 = apply(data_origscale_norm[,group_5a], 1, mean, na.rm = T) - 
                 apply(data_origscale_norm[,group_8a], 1, mean, na.rm = T))
# The values that are NA had no expression in any sample, so count as 0
D$D1 = na0(D$D1)
D$D2 = na0(D$D2)
# Average
D$Dav = sapply(1:nrow(D), function(i) mean(c(D$D1[i], D$D2[i])))                 
D = join(metabolite_annotations, D, type = 'right')
D = D[order(D$Dav, decreasing = T),]

#write.csv(D, 'Results/Tables/Calculation_natural_degradation.csv',
#            row.names = F, na = '')

pdf('Results/Plots/Histogram_D1_degradation.pdf',
    height = 5, width = 5)
hist(D$D1, 
     col = 'darkgrey',
     main = 'Degradation with stretch (D1)',
     las = 1,
     xlab = 'Degradation (D1)',
     ylab = 'Number of metabolites')
dev.off()
pdf('Results/Plots/Histogram_D2_degradation.pdf',
    height = 5, width = 5)
hist(D$D2, 
     col = 'darkgrey',
     main = 'Degradation without stretch (D2)',
     las = 1,
     xlab = 'Degradation (D2)',
     ylab = 'Number of metabolites')
dev.off()
pdf('Results/Plots/Histogram_average_degradation.pdf',
    height = 5, width = 5)
hist(D$Dav, 
     col = 'darkgrey',
     main = 'Average degradation (Dav)',
     las = 1,
     xlab = 'Degradation (Dav)',
     ylab = 'Number of metabolites')
dev.off()


# Calculations combined
calcs_combined = data.frame(Metabolon_Metabolite_ID = row.names(data_origscale_norm),
                            stringsAsFactors = F)
all_groups = unique(sample_annotations$Group_Number)
for (group_name in all_groups[!all_groups %in% c('11', '12')]){
  calcs_combined[paste0('group_', group_name, '_avg')] = get(paste0('group_', group_name, '_avg'))
}
calcs_combined = join(metabolite_annotations, calcs_combined)
calcs_combined = calcs_combined[order(calcs_combined$Sub_Pathway),]
calcs_combined = calcs_combined[order(calcs_combined$Super_Pathway),]

#write.csv(calcs_combined, 'Results/Tables/Average_metabolite_abundance_by_group.csv',
#          row.names = F, na = '')


# HEATMAP
row.names(calcs_combined) = calcs_combined$Biochemical_Name
heatmap_data = as.matrix(calcs_combined[,9:ncol(calcs_combined)])
heatmap_data[is.na(heatmap_data)] = 0
heatmap_data = heatmap_data[apply(heatmap_data, 1, function(x) !all(x==0)),]
colnames(heatmap_data) = gsub('group_|_avg', '', colnames(heatmap_data))

heatmap_palette = colorRampPalette(c('white', PLOT_COLS[11], PLOT_COLS[11]))(n = 100)
group_annotations = unique(sample_annotations[,c('Group_Number', 'Stretch',
                                                 'Donor_Type', 'Compartment',
                                                 'Condition')])
row.names(group_annotations) = group_annotations$Group_Number

stretch_colors = c('black', 'white')
names(stretch_colors) = c('Plus', 'Minus')
donor_colors = PLOT_COLS[c(8, 10)]
names(donor_colors) = c('EED', 'Healthy')
compartment_colors = c(PLOT_COLS[c(1, 3)], 'grey')
names(compartment_colors) = c('Basal', 'Apical', 'Cell Lysate')
condition_colors = PLOT_COLS[c(2,4,5,6,7)]
names(condition_colors) = c('Outflow', 'Inlet', 'Dry pellet',
                            'Loading media (DM)', 'Loading media (HBSS)')

calcs_combined2 = calcs_combined[calcs_combined$Biochemical_Name %in% row.names(heatmap_data),]
pathway_colors = PLOT_COLS[c(2,4,5,6,7,8,9,10)]
names(pathway_colors) = unique(calcs_combined2$Super_Pathway)
row_colors = matrix(nrow = 1, ncol = nrow(heatmap_data),
                    pathway_colors[calcs_combined2[row.names(heatmap_data), 'Super_Pathway']]
)

heatmap_data_amir_ordering = heatmap_data[,c('6b', '2b',  # Basal HBSS Healthy EED
                                             '10b',       # Basal Outflow Healthy -
                                             '4b',        # Basal Outflow EED -
                                             '5a', '1a',  # Apical DM Healthy EED
                                             '10a',
                                             '4a')] 

pdf(paste0('Results/Plots/Heatmap_average_metabolite_abundance_by_group_nocluster_subset.081821.pdf'),
    height = 15, width = 8)
heatmap.3(heatmap_data_amir_ordering[metabolite_subset,],
          main = paste0('Metabolite abundance across sample groups\n(Subset of ',
                        length(metabolite_subset), ' metabolites for which basal outflow > basal inlet)'),
          Rowv = T, 
          Colv = F,  
          ColSideColors = col_colors2,
          RowSideColors = row_colors2,
          keysize = 1,
          KeyValueName = 'Expression (scaled)',
          dendrogram = 'row',
          col = heatmap_palette,
          margins = c(3,10))
legend('topright',
       legend = c(unique(group_annotations$Donor_Type),
                  unique(group_annotations$Condition),
                  unique(group_annotations$Compartment)),
       pch = 15,
       col = c(donor_colors[unique(group_annotations$Donor_Type)],
               condition_colors[unique(group_annotations$Condition)],
               compartment_colors[unique(group_annotations$Compartment)]))
legend('topright',
       legend = names(pathway_colors),
       pch = 15,
       col = pathway_colors)
dev.off()
