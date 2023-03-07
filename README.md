# Circular Structural Model
The following includes all the MATLAB code used in my Master's Thesis for implementing maximum likelihood estimation of the circular structural model with the Von-Mises distribution 

## List of files

1. GetMoments - Computes the sample average of $X$ and $Y$, and $\widehat{\Sigma}$ the empirical covariance matrix of the pair $(X,Y)$.
2. UniformityTest - Function that apply the uniformity test and returns the $\chi^2(2)$ statistic and the corresponding p-values of the test.
3. CircleMLE_Uniform - Maximum likelihood estimation under uniformity of the observations across the circle.
4. FindMu - Estimating the 2 possible values of \$mu$ - the mean direction of the Von-Mises distribution using the spectal decomposotion of $\widehat\Sigma$.
5. bisection1 - Computing $\varkappa_L$, the lower bound of $\hat\varkappa$ using the bisection method.
6. get_mu_bounds1 - Lower and upper bounds for the 2 possible values of $\mu$.
7. VM_negLogLike3 - Function that computes the negative log-likelihood function of the sample $(x_i,y_i)$, $i=1,\dots,n$.
