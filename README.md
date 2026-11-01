# zilc
Zero-inflated linear-circular regression

R codes for the paper "Spline regression for zero-inflated circular data".

The files sim_1.R and sim_2.R contain the codes for the simulation studies, based on the ZIvM 1/ZIJP 1 models and the ZIvM 2/ZIJP 2 models, respectively. 

The files zivm1.R and vmsp.R contain the codes for the data analysis using the ZIvM models and the von Mises Spline (vMSp) model, respectively.

The comparison.R file contains codes for the model comparison with the Fisher-Lee, SPML, CIPC and GCPC models in terms of the circular CRPS and Leave-one-out Cross Validation.

Due to confidentiality restrictions, we do not have permission to make the original data analysed in our manuscript publicly available. As an alternative, a similar dataset (dat.csv) generated using parametric bootstrapping is provided here solely for illustrative purposes and for use with the accompanying code. The results obtained from analysing this dataset will naturally differ from those obtained using the original data.
