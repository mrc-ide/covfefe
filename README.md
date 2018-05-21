# covfefe
[![Travis-CI Build Status](https://travis-ci.org/mrc-ide/covfefe.svg?branch=master)](https://travis-ci.org/mrc-ide/covfefe)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/mrc-ide/covfefe?branch=master&svg=true)](https://ci.appveyor.com/project/mrc-ide/covfefe)
[![Coverage Status](https://img.shields.io/codecov/c/github/mrc-ide/covfefe/master.svg)](https://codecov.io/github/mrc-ide/covfefe?branch=master)
[![Documentation](https://img.shields.io/badge/documentation-click%20here!-blue.svg)](https://mrc-ide.github.io/covfefe/)

The R package *covfefe* performs Computational Validation For Estimating Falciparum Epidemiology (ahem). More specifically, it provides a flexible framework for simulating genetic data from *P. falciparum* transmission models.

We start with data on infection histories - i.e. who became infected when, who did this infection originate from, when did this episode recover, etc. We then specify a subset of samples that we are interested in - for example, 100 individuals at day 365. *covfefe* then prunes the infection history data down to a reduced set that contains only the events that directly impact on the samples that we are interested in. Often the number of sampled hosts is far smaller than the total infected population, or perhaps sampling only occurs at a few time points, meaning the pruned infection history can be *much* smaller than the complete infection history. Finally, *covfefe* simulates genetic data from the pruned infection history and from a set of user-defined parameters that specify the genetic model, for example the recombination rate and the number of loci.

By using the infection history as the starting point we ensure that multiple epidemiological models can be plugged into the same tool. Any malaria transmission model that generates infection history data in the format required by *covfefe* can be used to simulate genetic data, and this format has been designed to be flexible enough to encompass a wide range of models. As a starting point, *covfefe* contains a simple individual-based model that simulates malaria transmission in multiple demes (partially isolated subpopulations), with or without migration between demes.

**This package is still in development**, meaning methods are likely to change in the near future. That being said, here are instructions on how to install and run *covfefe* as it currently stands.

## Installation

*covfefe* relies on the Rcpp package, which requires the following OS-specific steps:

* Windows
    - Download and install the appropriate version of [Rtools](https://cran.rstudio.com/bin/windows/Rtools/) for your version of R. On installation, ensure you check the box to arrange your system PATH as recommended by Rtools
* Mac OS X
    - Download and install [XCode](http://itunes.apple.com/us/app/xcode/id497799835?mt=12)
    - Within XCode go to Preferences : Downloads and install the Command Line Tools
* Linux (Debian/Ubuntu)
    - Install the core software development utilities required for R package development as well as LaTeX by executing

            ```
            sudo apt-get install r-base-dev texlive-full
            ```

Next, ensure that you have devtools installed by running
```r
install.packages("devtools")
```
Finally, install the *covfefe* package directly from GitHub by running
```r
devtools::install_github("mrc-ide/covfefe")
library(covfefe)
```

## Example analysis

We start by creating a covfefe project. This project will hold all inputs and outputs used by the package, and is a convenient way of keeping everything we need in one place.
```r
myproj <- covfefe_project()
```
Running `names(myproj)` we can see all the elements that make up a *covfefe* project, including parameters, distributions, etc. We will talk about each of these elements later on.

Next we need some data on infection histories. For convenience, we will use the in-built individual-based model with default settings.
```r
mysim <- sim_indiv()
```
Again, running `names(mysim)` we can see the individual elements of `mysim`, including `H` - the human population size in each deme, `daily_counts` - the number of susceptible vs. singly-infected vs. doubly-infected (and so on) individuals on each day, and `infection_history` - the main output we are interest in. For those who want to understand the infection history format, details can be found *here* (**TODO - link to format**).

We continue by loading the infection history into our project.
```r
myproj <- assign_infection_history(myproj, mysim$infection_history)
```
Next, we need to specify how many human hosts we want to sample from.  We do this using a simple matrix or data frame with three columns: column 1 gives the time of sampling, column 2 gives the deme that we are sampling from, and column 3 gives the number of hosts (hosts are drawn at random from the pool of blood stage individuals at that point in time). Note that in column 3, a value of -1 indicates that all blood stage hosts should be sampled. For this example we will sample at day 100 and day 365, in a single deme (we only have one deme), and drawing all hosts at the first time point and 100 hosts at the second time point.
```r
samp <- data.frame(times = c(100, 365), demes = 1, hosts = c(-1, 100))
```
When we come to load this data into our project there are a few other arguments to be aware of. First, we must specify
```r
myproj <- sample_hosts(myproj, samp, demes, recover_pop_counts = FALSE, max_infections = 5)
```


