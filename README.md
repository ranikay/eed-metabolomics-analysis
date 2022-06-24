# Nutritional deficiency in an intestine-on-a-chip recapitulates injury hallmarks associated with environmental enteric dysfunction

**Authors:**

Amir Bein<sup>1,+,&</sup>, Cicely W. Fadel<sup>1,2,&</sup>, Ben Swenor<sup>1</sup>, Wuji Cao<sup>1</sup>, Rani K. Powers<sup>1,++</sup>, Diogo M. Camacho<sup>1,+++</sup>, Arash Naziripour<sup>1</sup>, Andrew Parsons<sup>1</sup>, Nina LoGrande<sup>1</sup>, Sanjay Sharma<sup>1</sup>, Seongmin Kim<sup>1</sup>, Sasan Jalili-Firoozinezhad<sup>1,3</sup>, Jennifer Grant<sup>1</sup>, David T. Breault<sup>2,4,5</sup>, Junaid Iqbal<sup>6</sup>, Asad Ali<sup>6</sup>, Lee A Denson<sup>7,8</sup>, Sean R. Moore<sup>9</sup>, Rachelle Prantil-Baun<sup>1</sup>, Girija Goyal<sup>1</sup>, and Donald E. Ingber<sup>1,10,11*</sup>

**Affiliations:**

<sup>1</sup> Wyss Institute for Biologically Inspired Engineering, Harvard University, Boston, MA 02115, USA

<sup>2</sup> Department of Pediatrics, Harvard Medical School, Boston, MA, 02115, USA

<sup>3</sup> Department of Bioengineering and iBB - Institute for Bioengineering and Biosciences, Instituto Superior Técnico, Universidade de Lisboa, Lisboa, 1049-001, Portugal

<sup>4</sup> Division of Endocrinology, Boston Children's Hospital, Boston, MA 02115, USA

<sup>5</sup> Harvard Stem Cell Institute, Harvard University, Boston, MA 02139, USA

<sup>6</sup> Department of Paediatrics and Child Health, The Aga Khan University, Karachi 74800, Pakistan

<sup>7</sup> Division of Gastroenterology, Hepatology, and Nutrition, Cincinnati Children's Hospital Medical Center, Cincinnati, OH, 45229-3026, USA

<sup>8</sup> Department of Pediatrics, University of Cincinnati College of Medicine, Cincinnati, OH, 45229-3026, USA

<sup>9</sup> Division of Pediatric Gastroenterology, Hepatology, and Nutrition, Department of Pediatrics, University of Virginia, Charlottesville,22908-0386, USA

<sup>10</sup> Harvard John A. Paulson School of Engineering and Applied Sciences, Harvard University, Cambridge, MA 02138, USA

<sup>11</sup> Vascular Biology Program and Department of Surgery, Harvard Medical School and Boston Children’s Hospital, Boston, MA 02115, USA


<sup>&</sup> These authors contributed equally

<sup>+</sup> Current address: Quris Technologies. Boston, MA, 02118, USA

<sup>++</sup> Current address: Pluto Biosciences, Inc., Golden, CO, 80402, USA

<sup>+++</sup> Current address: Rheos Medicines, Cambridge, MA, 02142, USA

<sup>*</sup> Corresponding author, `don.ingber@wyss.harvard.edu`


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
