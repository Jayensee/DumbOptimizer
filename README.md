# DumbOptimizer
Apps Script implementation of a fairly generic (and simple) stochastic optimizer and a weighting scheme that is robust to non-uniform sample densities/rates. Meant to be used by people who are comfortable working with Google Apps Scripts in Google Sheets.

## How it works
The optimizer starts with a user-defined initial estimation of the parameters to optimize and at every iteration tries random values around the best parameters it has found. Every time it finds an improvement it increases the reach of the randomization step, and decreases it when it doesn't. Finishes after a set number of iterations or when the reach of the randomization step gets low enough (whichever happens first).

## How to use
To use copy the code to an Apps Script file in the editor (Extensions > Apps Script on Google Sheets).
The documentation is in the code itself.

## What is in this?
The code contains 3 functions:
1. **dumb_optim**: stochastic optimizer
2. **compute_weights**: calculates kernel density with a Gaussian kernel then takes the reciprocals and normalizes them
3. **ldx_dx_model_error**: the only function to be optimized currently implemented 



### Important Note: whatever function you want to optimize you will need to implement and pass it onto the optimizer as an argument (unless by chance you happen to want the exact same function I wanted when I made this).

## Why is this repo even a thing?

The main reason I'm making this public is it took me forever to find a way to calculate the error in a way that isn't biased to points that were sampled close together.   
That and it's useful to have a generic optimizer that doesn't have a ton of random restrictions, even if it costs in efficiency (especially because apparently there isn't a single decent solver add-on for Google Sheets). 

If you need a way to fit whatever weird model you have with data from a spreadsheet this will probably get good enough results.    
If you took samples at irregular time intervals and every regression you try overfits to clusters of points check out the compute_weights function.

If you for some reason want to fit the parameters of a variable that decays exponentially but whose derivative increases as a latent variable decays, also exponentially, then that's also here (i.e. a chemical reaction where the product immediately starts reacting with something else).
