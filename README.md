# Nutritional deficiency recapitulates intestinal injury associated with environmental enteric dysfunction in patient-derived Organ Chips

Amir Bein^1,13^, Cicely W. Fadel^1,2,13^, Ben Swenor^1^, Wuji Cao^1^, Rani K. Powers^1^, Diogo M. Camacho^1,3^, Arash Naziripour^1^, Andrew Parsons^1^, Nina LoGrande^1^, Sanjay Sharma^1^, Seongmin Kim^1^, Sasan Jalili-Firoozinezhad^1,4^, Jennifer Grant^1^, David T. Breault^2,5,6^, Junaid Iqbal^7^, Asad Ali^7^, Lee A Denson^8,9^, Sean R. Moore^10^, Rachelle Prantil-Baun^1^, Girija Goyal^1^, and Donald E. Ingber^1,11,12*^

### Analysis of metabolomics data from patient-derived organ chips to investigate environmental enteric dysfunction

Install required R packages:

```
plyr
RColorBrewer
beeswarm
```


Run the scripts in order to perform analysis:

**1-compute_total_ion_count.R** - Calculates total ion count in outflow to determine whether it's a reasonable proxy for sample-wise variance (in lieu of Bradford for all samples). Normalizes each sample (origscale data) by its total ion count.

**2-calculate_metrics.R** - Calculates natural degradation, cellular utilization, and uptake and transfer. Saves the heatmap for Figure 4b.