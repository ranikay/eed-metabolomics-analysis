# --------------------------------------------------------------------------- #
# Script 1: Compute total ion count
#
# Part A: Calculate total ion count in outflow to determine whether it's a
# reasonable proxy for sample-wise variance (in lieu of Bradford for all samples)
#
# Part B: Normalize each sample (origscale data) by its total ion count and fill
# missing values with the minimum value for a metabolite
# --------------------------------------------------------------------------- #

library(plyr)
library(RColorBrewer)

### PART A

# Read sample data
sample_annotations = read.csv('Data/All_sample_data.csv', stringsAsFactors = F)

# Use just the outflow samples for this analysis
outflow_sample_annotations = sample_annotations[sample_annotations$Condition == 'Outflow',]
outflow_samples = outflow_sample_annotations$Sample_Name
row.names(outflow_sample_annotations) = outflow_samples

# Read metabolite data (origscale from Metabolon)
data_origscale = read.csv('Data/Media_metabolite_data_origscale.csv', 
                          stringsAsFactors = F, row.names = 1)
data_origscale_outflow = data_origscale[,outflow_samples]

# Calculate total ion count per sample
tic_data = data.frame(Sample_Name = names(data_origscale),
                      Total_Ion_Count_OrigScale = apply(data_origscale, 2, sum, na.rm = T))
tic_data = join(tic_data, sample_annotations)
tic_data = tic_data[order(tic_data$Total_Ion_Count_OrigScale, decreasing = F),]

write.csv(tic_data, 'Results/Tables/Media_total_ion_count_by_sample_all.csv', row.names = F, na = '')

### PART B

# Normalize each sample (origscale data) by its total ion count
row.names(tic_data) = tic_data$Sample_Name
data_origscale_norm = data_origscale
for (sample_name in names(data_origscale_norm)){
  data_origscale_norm[,sample_name] = data_origscale_norm[,sample_name] / tic_data[sample_name, 'Total_Ion_Count_OrigScale']
}

# Rescale each metabolite values by dividing each metabolite by its root mean square
data_origscale_norm = as.data.frame(cbind(Metabolon_Metabolite_ID = row.names(data_origscale_norm),
                                          data_origscale_norm))
data_origscale_norm2 = as.data.frame(t(scale(t(data_origscale_norm[,-1]), center = F)))
data_origscale_norm2 = as.data.frame(cbind(Metabolon_Metabolite_ID = row.names(data_origscale_norm2),
                                           data_origscale_norm2))

write.csv(data_origscale_norm2, 
          'Data/All_metabolite_data_origscale_tic_normalized_rescaled.csv',
          row.names = F, na = '')
