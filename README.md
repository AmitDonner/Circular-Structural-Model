# Circular Structural Model
The following includes all the MATLAB code used in my Master's Thesis for implementing maximum likelihood estimation of the circular structural model with the Von-Mises distribution 

## List of files

GetMoments - Computes the moment estimators $\bar{X}_n = \frac{1}{n}\sum_{i=1}^{n}X_i, \; \bar{Y}_n = \frac{1}{n}\sum_{i=1}^{n}Y_i$ and $\widehat{\Sigma}$ the empirical covariance matrix of the pair $(X,Y)$.
UniformityTest - Function that apply the uniformity test and returns the $\chi^2(2)$ statistic and the corresponding p-values of the test.
CircleMLE_Uniform - Maximum likelihood estimation under uniformity of the observations across the circle.
FindMu - Estimating the 2 possible values of \$mu$ - the mean direction of the Von-Mises distribution using the spectal decomposotion of $\widehat\Sigma$.
bisection1 - Computing $\varkappa_L$, the lower bound of $\hat\varkappa$ using the bisection method.
get_mu_bounds1 - Lower and upper bounds for the 2 possible values of $\mu$.
VM_negLogLike3 - Function that computes the negative log-likelihood function of the sample $\{(x_i,y_i),\;i=1,\dots,n\}$.
