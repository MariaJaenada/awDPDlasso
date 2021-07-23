# awDPDlasso
Code for fitting Linear and Logistic Regression Models, using weighted adaptative lasso, SCAD and MCP penalties based on the DPD loss procedure.

This code was used to obtain the regression estimators in the paper "Robust adaptive variable selection in ultra-high dimensional linear regression models" (arXiv:2004.05470) and 
Robust adaptive Lasso in high-dimensional logistic regression with an application to genomic classification of cancer patients.

The files hgic_dpd_linear and hgic_dpd_logit contains the main functions, awDPDlasso.lm and awDPDlasso.logistic to fit the linear and the logistic model respectively.
In awDPDlasso.lm, the user should introduce as arguments the data (matrices X and Y), an initializer mode (between RLARS or standar DPD-LASSO), a weight.mode (between "lasso", "adaptive" and "scad" for the three differents weights functions used in the main paper), a penalty.mode (between "lasso", for lasso and adaptive methods, "scad" or "mcp") and a mode of choosing the lambda grid. 
In awDPDlasso.logistic, the user should introduce as arguments the data (matrices X and Y), a weight.mode ("between "lasso", "adaptive" and "scad") and a mode of solving the auxiliary linear regression problem ("LS" or "DPD").

In both cases, "lambda0" uses an auxiliar function to determine the extremes of the path, and "given" uses the "lmin" and "lmax" limits provided by the user. Finally, "nlambda" determines the length of the lambda grid. The tunning parameter (lambda) selection is performed using the HBIC criterion (see the main paper for more details).

This code was partially based on the code implemented in https://github.com/cran/gamreg for Robust regression via gamma-divergence with L1, elastic net and ridge.
