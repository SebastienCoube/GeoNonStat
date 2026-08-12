We use the funtion plotPointillistPainting() in order to visualize the response variable.

```{r}
plotPointillistPainting(dataset$locs, dataset$observed_field, cex = 8, main = "response")
```

Not very conclusive. How about one variable from X?
  
  ```{r}
plotPointillistPainting(dataset$locs, dataset$X[,1], cex = 8, main = "X[,1]")
```
We can see that the first variable from X is equal to the first spatial coordinate. 

## Initializing the model

### Vecchia approximation

The Vecchia approxiation is used to build the Gaussian Processes used in the model. 

```{r}
vecchia_approx <- createVecchia(dataset$locs, m=6)
```
### Predictive Processes (PP) 
Predictive Processes are used to model spatial variation in nonstationary model parameter. Aside of the Vecchia approximation, a positive spatial range parameter is needed. 
Let's make one PP for the nonstationary model range, and another for the nonstationary model noise. 


```{r}
PP_range = createPP(vecchia_approx, matern_range = .5)
```
```{r}
PP_range = createPP(vecchia_approx, matern_range = .2)
```

### Model initialization with GeoNonStat function

The model is intialized using the GeoNonStat function.
```{r}
geo_non_stat = GeoNonStat(
  observed_field = dataset$observed_field, # response
  vecchia_approx = vecchia_approx, # Vecchia approximation
  X = dataset$X, # explanatory variables for the response 
  noise_X = dataset$noise_X, noise_PP = PP_noise,  # explanatory variables and PP for the nonstationary noise variance 
  range_PP = PP_range # PP for the nonstationary range 
)
```
The resulting object has its own class and print function. 
```{r}
class(geo_non_stat)
print(geo_non_stat)
```
## Running the model 

The model is run using the future setup. Each chain can have several threads.
```{r}
 run_time = Sys.time()
 future::plan(strategy = "multisession", workers = 2)
 geo_non_stat <- GeoNonStatMcmc(
   geo_non_stat, n_iterations = 200, 
   n_chains_in_parallel = 2, n_threads_per_chain = 3)
 run_time = run_time - Sys.time()
```

```{r}
mcmcDiags(geo_non_stat, burn_in = .5)
```
```{r}
tracePlots(geo_non_stat, burn_in = .5)
```








































