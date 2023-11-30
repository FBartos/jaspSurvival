Non-parametric Survival Analysis
==========================

The non-parameteric survival analysis allows the user to estimate the survival probabilities and test the null hypothesis that the survival curves are equal. 

### Input
-------

#### Assignment Box 
- Time to Event: In this box the survival time variable is selected.  
- Event Status: In this box the variable containing the event indicator is selected. 
- Event Indicator: This dropdown defines which level of the Event Status variable corresponds to an event.
- Factors: In this box the factors for stratifying the analysis are selected.

#### Tests 
- Log-rank (Mantel-Haenszel): Log-rank (Mantel-Haenszel) test 
- Peto and Peto: Peto and Peto's test 
- Flemming-Harrington: Flimming-Harrington test

#### Life table
Summarizes survival using a life table.

#### Survival curve plot
Displays survival curves plot.
- Confidence interval: Whether a 95% CI should be visualized in the survival curve plot.
- Risk Table: Displays table summarizing the number at risk at a given time.
- Cummulative events tale: Displays table summarizing the cumulative number of events at a given time.
- Censoring plot: Displays distribution of the censored observations in each strata.
  - Cumulative: Changes the figure to a table summarizing the cummulative number of censored observations at a given time.
- Legend: Location of the legend.
- Color pallete: Color pallete for coloring strata.

### Output
-------

#### Kaplan-Meier Summary Table
- Strata: The strata for which the corresponding summary is generated. 
- N: Total number of observations. 
- Events: Number of events. 
- Restricted Mean: Restricted mean survival estimate.
- Standard Error: Standard error of the restricted mean estimate.
- Median Survival: Median survival.
- 95% CI for the restricted mean. 
    - Lower: The lower bound of the confidence interval. 
    - Upper: The upper bound of the confidence interval. 

#### Tests Table
- Test: Test's name. 
- Chi Square: The value of the Chi square test statistic. 
- df: Degrees of freedin if the Chi square test.
- p: The p-value. 


### References
-------
- Mills, M. (2010). Introducing Survival and Event History Analysis. SAGE.

### R-packages
---
- survival
- survminer

### Example 
--- 
- For an example go to `Open`--> `Data Library` --> `Surival` --> `Leukemia`. 


