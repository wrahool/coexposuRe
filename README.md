# coexposuRe
coexposuRe is an R package that helps simulate audience coexposure networks. It allows you to model audience behavior in an artificial media environment by tuning various parameters that mimic certain constraints or affordances of real-world media environments. This in turn provides novel possibilities of testing theories related to audience behavior using formal network theoretic tools.

The ideas underpinning `coexposuRe` are available in this paper:

> Mukerjee, S. (2021). A systematic comparison of community detection algorithms for measuring selective exposure in co-exposure networks. *Scientific Reports*, 11(1), 1-11.

To install `coexposuRe`, run `install_github("wrahool/coexposuRe")` in the R console.
The `install_github` function is available in the `devtools` package, which needs to be installed by running `install.packages("devtools")` and included using `library(devtools)`.

`coexposuRe` provides two main functions:

1. `simulate_single_network`: generates a single co-exposure network by specifying certain parameter values
2. `analyze_simulated_networks`: generates multiple networks with the same parameters as above, but with different values of the randomizing parameter, and then analyzes them. Currently, two analyses are supported, community detection and centralization.

For details type `?simulate_single_network` and `?analyze_simulated_networks`

## Citation

If you use `coexposuRe` for your work, please consider citing the following paper:

```
@article{mukerjee_systematic_2021,
	title = {A systematic comparison of community detection algorithms for measuring selective exposure in co-exposure networks},
	volume = {11},
	url = {https://www.nature.com/articles/s41598-021-94724-1},
	doi = {10.1038/s41598-021-94724-1},
	journal = {Scientific Reports},
	author = {Mukerjee, Subhayan},
	year = {2021},
}
```
