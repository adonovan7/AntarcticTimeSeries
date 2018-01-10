# Time Series Analysis of Mean Monthly Temperatures in Antartica 

## Overview

I looked at monthly surface temperature means from March of 1944 to November of 2012 at the Faraday station in Antarctica, with the ultimate goal of accurately forecasting future temperature fluctuations

To model the data, I used Box-Jenkins seasonal and non-seasonal autoregressive integrated moving average (SARIMA and ARIMA) models to capture the changes in the average monthly measurements, properly process any homogenous non-stationarity that might be present (small clusters of similarly-behaving data points), and forecast future temperature data points


## Conclusions

After a in-depth analysis of the time series data as well as evaluation of a large selection of potential representative models, I chose a SARIMA (8,0,1) x (0,1,1)12 model to fit the data.

The model did not pass the McLeod-Li test for squared residuals—however—the model did pass the rest of the diagnostic tests. In terms of forecasting, the model clearly produced forecasts that were slightly shifted up from where the actual data left off. This could either be due to a fault in the indexing or the model itself. Despite the offset, the forecast seems to follow the sinuous shape of the real data very well. 

On a different note, because there was no upward linear trend present after differencing for seasonality, I was not able to conclude that the mean monthly surface temperatures recorded at the station have increased significantly over the past 70 years. However, more recent data as well as data from other stations in the South Pole could be extremely useful in furthering such research. 


## Acknowledgments

This project was done for Professor Raya Feldman's PSTAT 174/274 Time Series class at UC Santa Barbara. I heavily used her lecture notes in addition to outside online and textual sources, including: 

*  Robert Shumway and David Stoffer’s textbook Time Series Analysis and its Applications, 3rd Edition 
* [O Texts](https://www.otexts.org/fpp/8/1)
* [Duke Lecture Notes](http://people.duke.edu/~rnau/411arim3.htm)
* [A Little Book of R for Time Series](https://a-little-book-of-r-for-time-series.readthedocs.io/en/latest/src/timeseries.html)

## Data Sources:

* [Original Source](www.bas.ac.uk)
* [Easy to Access Version](http://www.nerc-bas.ac.uk/icd/gjma/faraday.temps.html)