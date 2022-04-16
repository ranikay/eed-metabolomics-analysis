# Nutritional deficiency recapitulates intestinal injury associated with environmental enteric dysfunction in patient-derived Organ Chips

**Authors:**

Amir Bein<sup>1,13</sup>, Cicely W. Fadel<sup>1,2,13</sup>, Ben Swenor<sup>1</sup>, Wuji Cao<sup>1</sup>, Rani K. Powers<sup>1</sup>, Diogo M. Camacho<sup>1,3</sup>, Arash Naziripour<sup>1</sup>, Andrew Parsons<sup>1</sup>, Nina LoGrande<sup>1</sup>, Sanjay Sharma<sup>1</sup>, Seongmin Kim<sup>1</sup>, Sasan Jalili-Firoozinezhad<sup>1,4</sup>, Jennifer Grant<sup>1</sup>, David T. Breault<sup>2,5,6</sup>, Junaid Iqbal<sup>7</sup>, Asad Ali<sup>7</sup>, Lee A Denson<sup>8,9</sup>, Sean R. Moore<sup>10</sup>, Rachelle Prantil-Baun<sup>1</sup>, Girija Goyal<sup>1</sup>, and Donald E. Ingber<sup>1,11,12*</sup>

**Affiliations:**

<sup>1</sup> Wyss Institute for Biologically Inspired Engineering, Harvard University, Boston, MA 02115, USA
<sup>2</sup> Department of Pediatrics, Harvard Medical School, Boston, MA, 02115, USA
<sup>3</sup> Current address: Rheos Medicines, Cambridge, MA, USA
<sup>4</sup> Department of Bioengineering and iBB - Institute for Bioengineering and Biosciences, Instituto Superior Técnico, Universidade de Lisboa, Lisboa, Portugal
<sup>5</sup> Division of Endocrinology, Boston Children's Hospital, Boston, MA 02115, USA
<sup>6</sup> Harvard Stem Cell Institute, Harvard University, Boston, MA 02139, USA
<sup>7</sup> Department of Paediatrics and Child Health, The Aga Khan University, Karachi 74800, Pakistan
<sup>8</sup> Division of Gastroenterology, Hepatology, and Nutrition, Cincinnati Children's Hospital Medical Center, Cincinnati, OH, USA.
<sup>9</sup> Department of Pediatrics, University of Cincinnati College of Medicine, Cincinnati, OH, USA.
<sup>10</sup> Division of Pediatric Gastroenterology, Hepatology, and Nutrition, Department of Pediatrics, University of Virginia, Charlottesville, USA.
<sup>11</sup> Harvard John A. Paulson School of Engineering and Applied Sciences, Harvard University, Cambridge, MA 02138, USA.
<sup>12</sup> Vascular Biology Program and Department of Surgery, Harvard Medical School and Boston Children’s Hospital, Boston, MA 02115, USA.
<sup>13</sup> These authors contributed equally: Amir Bein, Cicely W. Fadel


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