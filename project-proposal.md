## Project Progress Report

### Theoretical Motivation

Diagnosing Autism Spectrum Disorder (ASD) is a very difficult task as there isn't a full proof method like a blood test or brain MRI Scan to diagnose the disorder. Doctors look at the child’s behavior and development to make a diagnosis.[1](https://www.cdc.gov/ncbddd/autism/screening.html). There are many children that are misdiagnosed. For example, it has been argued, overall better social and linguistic skills of ASD women have contributed to a misdiagnosis or late diagnosis. The aim of this project is to find biomarkers using the fMRI data of typically developing(TD) people and people with ASD that can help predict if a person suffers from ASD or not. We will be working with ABIDE I data, specifically the NYU's site resting state fMRI data


Moreover ASD women profile has not been fully investigated. Possibly, due to striking differences in prevalence (3:1) and misdiagnosis. We wish to look into that as well if time permits.

### Research Design
We plan to start without a hypothesis, with a completely data based approach.
1. Parcellate the whole brain using some good Atlas.
2. Calculating the pair-wise Functional Connectivity using average BOLD time series from all these regions got using Atlas.
3. Then doing a group level analysis to see the FC different in ASD vs TD.
4. Constructing a criterea for single subject diagnosis of ASD.
5. Doing analysis of the two groups of FC maps using interpretable Machine Learning approaches.
6. Comparing the approaches {3,4} vs {5}.

### Statistical Analysis

#### Standard statistical testing approach
1. Define ROI's (Regions of Interests).    
    - Atlas based  
        - Pros:
          - By far, I think this is the easiest approach.  
        - Cons:
          - Not recommended in longitudinal studies.  

  2. Extract time series of ROIs by taking average of all voxel time series.
    - Pros:
      - Easy to implement

    - Cons:
      - Might introduce artifacts.

  3. Run pearson and Partial correlations to get the FC matrices.  

#### Machine learning approach

There are three sequential steps:

   1. Apply feature selection algorithms to select subsets of voxels for use in model construction:
      - Multivariate feature selection algorithm
      - Recursive Feature Elimination (RFE):
        - Uses the training algorithm (support vector machine) recursively to eliminate irrelevant voxels and estimate informative spatial patterns. [Article](http://www.ncbi.nlm.nih.gov/pubmed/18672070)    

   2. Buid predictive models to classify the two classes i.e ASD men vs ASD women, then compare the accuracy of the models.

   3. Finding discriminating voxels and region

      - Voxels which lead to the discrimination between the classes are identified.
      - Finding the regions based on the discriminating voxles with best predictive accuracy.

At last,we can compare the results of statistical and Machine learning approach.

### Code Development
- We will be primarily working with Nipype.
- We have converted the data to BIDS format.
- Currently we are working on constructing a preprocessing pipeline.
