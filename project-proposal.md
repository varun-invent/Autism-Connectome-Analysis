## Project Proposal
Following is the project proposed by Laura V  [[Link to Github Repo](https://github.com/lau-v/Autism-Connectome-Analysis)]. I just created the .md version of it.

Hi everyone!
Following the advice of Rapid Rounds for Reproducible Science:
[Youtube Link](https://youtu.be/huZxlkAwRNw), I brainstormed some ideas for a suitable research project, based on the initial idea we used to pitch this OBI project. 
Please, feel free to provide feedback in any of these domains. This is a very rough idea in need of refinement. 

### Theoretical Motivation
ASD women profile has not been fully investigated. Possibly, due to striking differences in prevalence (3:1) and misdiagnosis. It has been argued, overall better social and linguistic skills of ASD women have contributed to a misdiagnosis or late diagnosis. The aim of this project is to test whether these skills indeed differ between men and women. [Currently, I am looking whether there is previous research investigating phenotypic differences considering these two variables].

### Research Design
To study social and linguistic networks of ASD a mean difference in the functional connectivity of the main social and linguistic centers will be addressed. Handedness, IQ, and motion artifacts (mean frame-wise displacement) will be used as confounding factors.  
- [] @Eve02 I am not sure what was your plan to incorporate differences in short range vs long range connectivity. Were you thinking of age differences as well? Maybe, you could provide some input here.  

To proceed with the analysis, I think of two suitable ways. Either:
1.  We contrast ASD groups with their neurotypical controls (i.e. ASD men vs neurotypical men, ASD women vs ASD neurotypical women).  
  -  Pros:  
      -Better overview of sex differences, useful for interpreting results.  
      -To really assume differences are ASD intrinsic, it is better to first test ASDs vs controls.   
  -  Cons:  
      -Less features will be available to compare between groups (only handedness, IQ and motion artifacts).   
       As far as I know, no other psychological batteries were applied to both (ASD vs controls).  
       
2.  We contrast solely group differences between ASD men and women (i.e. ASD men vs ASD women).   
  -  Pros:     
      -Characteristic contributing to these differences might be addressed by symptom severity only ASD population may have a                   record of: Ex. age of diagnosis, higher thresholds in sensory profiles, better communication and social skills as assessed               by scales used in the diagnostic instruments for ASD (ADOS and ADI).   
  -  Cons:      
      -Differences in group sizes will be huge. Perhaps, by matching subjects this could be somehow be solved?   
      -Lack of a baseline group (neurotypical controls)       
      -On a second thought, maybe it is better to start by contrasting between group differences rather than within group differences.  

### Statistical Analysis
[Since the original proposal includes a machine learning approach, I assume after standard pipeline preprocessing, the analysis requires two stages].

1.  Standard pipeline preprocessing.I think at this point, this would be the main issue. Decide how are we going to proceed for:
 [ [pdf with pipeline preprocessing using FSL] (http://www.humanbrainmapping.org/files/2015/Ed%20Materials/FSL_PreProcessing_Pipeline_OHBM15_Jenkinson.pdf)]
  - [] Distortion correction
  - [] Spatial smoothing
  - [] Temporal filtering
  - [] GLM
  - [] Resampling to standard space.  
  
2.  Statistical Analysis.
 -  Social and linguistic networks could be defined using as seed regions main social and linguistic centers (extracted from available   literature).  
 
  GLM approach  
  1. [] Define ROI's (Regions of Interests).  
      Two possible approaches:  
        - Atlas based  
          -Pros: By far, I think this is the easiest approach.  
          -Cons: Not recommended in longitudinal studies.  
         - ICA parcellation  
           -Pros: Preferred using subjects template.  
           -Cons: Does anyone have experience with this?  
  2. [] Extract time series of ROIs
  3. [] Run pearson and Partial correlations to test sex differences in ASD population.
 
### Machine learning approach
   
   There are three sequential steps that i can suggest
  
   3.1 Apply feature selection algorithms to select subsets of voxels for use in model construction ,one is multivariate             feature selection algorithm,Recursive Feature Elimination (RFE) that uses the training algorithm (support vector             machine) recursively to eliminate irrelevant voxels and estimate informative spatial patterns.
        `Article://www.ncbi.nlm.nih.gov/pubmed/18672070`    
   3.2 Buid predictive models to classify the two classes i.e ASD men vs ASD women, then compare the accuracy of the models. 
   
   3.3 Finding discriminating voxels and region
   
    3.3.1 Voxels which lead to the discrimination between the classes are identified.
    3.3.2 Finding the regions based on the discriminating voxles with best predictive accuracy.
    
`In the last,we can compare the results of statistical and Machine learning approach.`

### Code Development
I have not made any extensive research in this part yet, though we can base our code on available nypipe projects.

### Other Issues:
I noticed we have to register to download the database. Have any of you completed the registry? I wondered whether we can archive the data online. Otherwise it would occupied large memory storage on my computer and I guess it would be computationally intense (thus, not recommended).
