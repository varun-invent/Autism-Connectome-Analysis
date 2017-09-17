## Project Progress Report

### Theoretical Motivation

Autism Spectrum Disorder (ASD) is a neurodevelopmental disorder characterized by deficits in social communication and social interaction, and by restrictive repetitive patterns of behaviors. ASD heterogeneity is huge, symptom severity varies from severe impairment to mild impairment, and from requiring substantial support to requiring support. Consequently, a categorial approach and a behavioral assessment of the ASD core symptoms fail to detect ASD subjects in the upper spectrum. ASD subjects with mild impairments are susceptible of a misdiagnosis or a late diagnosis [1](https://www.ncbi.nlm.nih.gov/pubmed/25193140).

This project aims to find a connectivity biomarker that can use a resting state fMRI scan to differentiate a subject suffering from ASD from a subject with neurotypical development (TD). For that purpose, we will use resting state fMRI scans of ABIDE I NYU database. We aim to build a classifier with high sensitivity, but not at the expense of high specificity. However, we know it might be necessary to make a trade-off. 

### Research Design
We will employ a data driven approach. We plan to use Machine Learning to do the following:

1. Atlas based parcellation of the whole brain.
2. Extract the average BOLD time series of all those regions and calculate pairwise functional connectivity (FC) of those regions.
3. Run a group level analysis to contrast functional connectivity across groups (ASD vs. TD)
4. Use interpretable Machine Learning approaches to classify our groups according to its functional connectivity maps (ASD vs. TD).
5. Compare the aforementioned approaches in terms of sensitivity and specificity across regions.

#### Machine learning implementation

   To build our model, we will:
   1. Apply feature selection algorithms to select a subset of voxels using:
      - A multivariate feature selection algorithm, and 
      - A Recursive Feature Elimination (RFE):
        - We will use the training algorithm (support vector machine) recursively, to eliminate irrelevant voxels and estimate informative spatial patterns [2](http://www.ncbi.nlm.nih.gov/pubmed/18672070).    

   2. Buid predictive models to differentiate ASD subjects from TD subjects. 
   
   3. Compare the accuracy of the models above.

   4. Find discriminating voxels and regions

      - Discriminating voxels (those voxels with the highest predictive accuracy) will lead to class identification (ASD vs TD).
      - Discriminating regions will be found using discriminating voxels 
      
### Code Development
- We will be using mostly, Nipype.
- We have already converted our dataset to BIDS format.
- Currently, we are working in the preprocessing of our pipeline.

### References
1. Anderson, G. M. (2015). Autism biomarkers: challenges, pitfalls and possibilities. _Journal of autism and developmental disorders, 45(4)_, 1103-1113.

2. De Martino, F., Valente, G., Staeren, N., Ashburner, J., Goebel, R., & Formisano, E. (2008). Combining multivariate voxel selection and support vector machines for mapping and classification of fMRI spatial patterns. _Neuroimage_, 43(1), 44-58.

