### Future work (GWAEN)

#### Autocorrelated errors
The adaptive elastic net has been implemented and validated by simulation for Gaussian data with an iid error structure. One crucial feature of spatial data is that the errors are autocorrelated. To this point, simulations have indicated that the performance of the  GWAEN suffers when autocorrelated errors are introduced.

Typically, the "optimal" bandwidth is smaller when autocorrelated errors are introduced, perhaps because at small bandwidths, autocorrelated errors are difficult to distinguish from an intercept. This causes a reduction in the residual error at small bandwidths, compared to what it would be if the errors were iid. It is notable that the traditional GWR exhibits this same behavior.

Adapting the GWAEN to the setting of an autocorrelated error structure is one necessary enhancement.

#### Non-gaussian data
To date, validation of the GWAEN has been for gaussian data. The extension to any exponential-family distribution is straightforward: since the selection procedure for the GWAEN is based on the penalized log-likelihood, the same procedure should work as well in principle for any exponential family. Preliminary simulations bear out this judgement but the generalized GWAEN still lags the development of the gaussian GWAEN.

