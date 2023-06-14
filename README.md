# Circular Structural Model
The listed below files contains all the code needed for implementing the ML estimation of circles from noisy data, as described in the attached paper 

## List of files
1. GetMoments - Computes the sample average of $X$ and $Y$, and $\widehat{\Sigma}$ the empirical covariance matrix of the pair $(X,Y)$.
2. UniformityTest - Function that apply the uniformity test and returns the $\chi^2(2)$ statistic and the corresponding p-values of the test.
3. CircleMLE_Uniform - Maximum likelihood estimation under uniformity of the observations across the circle.
4. FindMu - Estimating the 2 possible values of \$mu$ - the mean direction of the Von-Mises distribution using the spectal decomposotion of $\widehat\Sigma$.
5. bisection1 - Computing $\varkappa_L$, the lower bound of $\hat\varkappa$ using the bisection method.
6. get_mu_bounds1 - Lower and upper bounds for the 2 possible values of $\mu$.
7. VM_negLogLike3 - Function that computes the negative log-likelihood function of the sample $(x_i,y_i)$, $i=1,\dots,n$.
8. CircleMLE - Wrapper of functions 1-7 that given a sample $(x_i,y_i)$, $i=1,\dots,n$ return an estimates of the circle parameters along with the nuisance paramteres.

## Instructions
Scripts 1-8 need to be downloaded to your working directory. Then, apply the function CircleMLE (i.e. script 8) which takes 2 arguments, $x$ values (vector) and $y$ (vector of the same length). 

CircleMLE returns a vector of length 6, with the following order - kappa, mu, rho, a, b, sigma. This function is a wrapper which uses function 1-7 in order to estimate the vector of parameters. 

There is no need to use functions 1-7 seperately. 

### Examples

* sol=CircleMLE(x,y) - Returns a vector of length 6, with the estimated parameters, ordered by kappa, mu, rho, a, b, sigma.
* [sol, P_val] = CircleMLE(x,y) - *sol* is the vector of estimated parameters and *P_val* is the P-value of the uniformity test.
* [sol, P_val, stat] = CircleMLE(x,y) - sol is the vector of estimated parameters, and *P_val* and *stat* are the uniformity test P-value and $\chi^2$ statistic respectively.


For any issue - feel free to contact me at *adonner@campus.haifa.ac.il*

