## Project Proposal
Following is the project proposed by Laura V  [[Link to Github Repo](https://github.com/lau-v/Autism-Connectome-Analysis)]. I just created the .md version of it.


Hi everyone!
Following the advice of Rapid Rounds for Reproducible Science:
[Youtube Link](https://youtu.be/huZxlkAwRNw), I brainstormed some ideas for a research project, based on the sample and approach initially proposed to pitch the idea for this OBI project. 
Please, feel free to provide feedback in any of these domains. It is more than welcome and needed to refined this very rough idea.


### Theoretical Motivation
Phenotypic traits of ASD women have not been fully addressed. Possibly, due to different prevalent rates (3:1) and misdiagnosis of ASD women. It has been argued, overall better social and linguistic skills of ASD women have contributed to its late diagnosis and under sampling. The aim of this project is to test whether these skills indeed differ between men and women. [Currently, I am looking whether there is previous research investigating phenotypic differences considering these two variables]

### Research Design
To study social and linguistic networks of ASD a mean difference in the functional connectivity of the main social and linguistic centers will be addressed. Handedness, IQ, and motion artifacts (mean frame-wise displacement) will be used as confounding factors.

Here, I think of two suitable ways. Either:
1. We contrast solely group differences between ASD men and women
(i.e. ASD men vs ASD women)
Pros:  
  - Characteristic contributing to these differences might be addressed by symptom severity only ASD population may have a record of: age of diagnosis, higher thresholds in sensory profiles, better communication and social skills as assessed by scales used in the diagnostic instruments for ASD (ADOS and ADI).
Cons:
  - Differences in group sizes will be huge. Perhaps, by matching subjects this could be somehow be solved?
  - Lack of a baseline group (neurotypical controls)

2. We contrast ASD groups with their neurotypical controls:
(i.e. ASD men vs neurotypical men, ASD women vs ASD neurotypical women)  
Pros:  
  - Better overview of sex differences, useful for interpreting results
  - To really assume differences are ASD intrinsic it is better to first test ASDs vs controls

  Cons:
  - Less features will be available to compare between groups (only handedness, IQ and motion artifacts). As far as I know, no other psychological batteries were applied to both (ASD vs controls). 
  
### Statistical Analysis
[Since the original proposal includes a machine learning approach, I assume the analysis requires two stages]
  -Social and linguistic networks could be defined using as seed regions main social and linguistic centers (extracted from available   literature).
  
1. Standard Analysis  
  1.1 Standard pipeline preprocessing 
  1.2.1 Define ROI (Regions of Interests):
      -Atlas based. 
      Pros: By far, I think this is the easiest approach. 
      Cons: Not recommended in longitudinal studies
      -ICA parcellation. Does anyone have experience with this?
  1.2.2 Pearson and Partial correlations will be used to test sex differences in ASD population.
    1.2.2 Extract time series of regions of interest (ROI), and run the corresponding analysis
 

2. Machine learning approach  
  - Here, input is needed on how to proceed with the approach of machine learning

### Code Development
I have not made any extensive research in this part yet, though we can base our code on available nypipe projects.

### Other Issues:
I noticed we have to register to download the database. Have any of you completed the registry? I wondered whether we can archive the data online. Otherwise it would occupied large memory storage on my computer and I guess it would be computationally intense (thus, not recommended).
