# Dynamic modeling via autoregressive conditional GB2 for cross-sectional maxima of financial time series data

This README provides demo codes for analyzing S&P 100 data as presented in the paper:

Fan, N., Zhang, C.M., and Zhang, Z.J. (2024). Dynamic modeling via autoregressive conditional GB2 for cross-sectional maxima of financial time series data. Journal of Business & Economic Statistics. [https://doi.org/10.1080/07350015.2025.2450488] 


Data: Stored in Microsoft Excel documents, in CSV format.

Tables:

•	Table 3: Generated using the MATLAB code ‘maximaplots_sp100.m’. 

•	Table 4 and Table 6: 
1.	Run the R codes ‘SP100-processed2022acf.R’ and ‘SP100-processed2022acgb2.R’ to obtain estimated parameter values, saved as:
"SP100paraacgb2New.csv”, 
"SP100likeacgb2New.csv"
2.	Use the MATLAB codes ‘Example02solveAcGB2SP100Tab6exp.m’ and
‘Example02solveAcFSP100Tab6exp.m’ to refine the estimates and compute Fisher information matrices and standard errors (s.e.). 

•	Table 8: Run the MATLAB code ‘Example02threedatasetsfittingKSgof.m’. 
	
  - Note: This performs Monte Carlo tests. The output may vary slightly from the values in the paper due to differences in seed numbers or random number generators.

•	Table 9: Run the MATLAB code ‘recoveredplots_sp100.m’.

•	Table 10: Run the MATLAB code ‘Example02threedatasetsforecasting.m’. 

  - Note: Same considerations as for Table 8 regarding output variability. 


Figures:

•	Figure 3: Run the MATLAB code ‘maximaplots_sp100.m’.

•	Figure 5: Run the MATLAB code ‘recoveredplots_sp100.m’.

•	Figure 8: Run the MATLAB code ‘Example02threedatasetsforecasting.m’.

•	Figure 11: Run the MATLAB code ‘recoveredplots_sp100.m’.

