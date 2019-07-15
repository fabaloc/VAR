# VAR Analyis

The Matlab code 'VARanalysis.m' run a Vector Autoregressive Analysis of US quarterly data from 1948Q1 to 2015Q4 : real gdp in QoQ growth rates (gdp),
QoQ inflation (inf) and the three months Treasury Bill rate (str).

VAR models with lags from 1 to 8 are estimated and then AIC and BIC are computed to find the optimal amount of lags to use.
Subsequently, a VAR(2) model is estimated and the stationarity of the process and absence of correlation among the residuals are established.


