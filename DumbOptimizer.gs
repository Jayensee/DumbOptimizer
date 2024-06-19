/**
 * Calculates the normalized reciprocal of the kernel density with a Gaussian kernel. The output weights compensate for non-uniform sample density. (Points close together get lower weight than isolated points so that the weight density is roughly constant).
 *
 * @param {array[n]} x The list of values along the axis in which constant density is desired.
 * @param {number} bandwidth The bandwidth for the Gaussian kernel.
 * @return {array[n]} weights The weights equal to the normalized reciprocal of the kernel density. 
 * @customfunction
 */
function compute_weights(x=[1, 1.1, 1.2, 4, 4.1, 4.2, 7], bandwidth=null) {
  // flatten x because sheets is annoying
  x = x.flat();

  // reasonable default bandwidth
  if (bandwidth == null) bandwidth = 2*(Math.max(...x)-Math.min(...x))/x.length;

  
  weights = new Array(x.length).fill(0);

  // calculate the density with a Gaussian kernel
  for (let i = 0; i < x.length; i++) {
    for (let j = 0; j < x.length; j++) {
      weights[i] += Math.exp(-0.5 * Math.pow((x[j] - x[i]) / bandwidth, 2));
    }
  }
  
  // set the weights as the reciprocal of the density (use the weights to compensate for the density)
  weight_sum = 0;
  for (let i = 0; i < x.length; i++) {
    weights[i] = 1/weights[i];
    weight_sum += weights[i];
  }

  // normalize the weights
  for (let i = 0; i < x.length; i++) {
    weights[i] = weights[i]/weight_sum;
  }
  
  return weights
}


/**
 * Calculates the error for a model where the output variable y is defined by the following system of differential equations: y'=-scale*u'-rd*(y-baseline), u'=-rl*u, y(0)=scale*id+baseline, u(0)=1. It's basically just the sum of two exponentials. 
 * Outputs the quadratic mean (also known as the root mean square) of the error weighted by the reciprocal of the kernel density across time.
 *
 * @param {array[5]} x The constants that define the model [rl, id, rd, scale, baseline].
 * @param {[array[n,3]]} fun_args The data points followed by the pre-computed weights.
 * @return {array[n]} error The quadratic mean of the error weighted by the reciprocal of the kernel density across time. 
 */
function ldx_dx_model_error(x, fun_args) {
  // flatten x because sheets is annoying
  x = x.flat();
  il = 1;
  rl = x[0];
  id = x[1];
  rd = x[2];
  scale = x[3];
  baseline = x[4];
  error = 0;

  weights = fun_args[1];

  for ([index, datum] of fun_args[0].entries()) {
    t = datum[0];
    y = datum[1];

    // calculate the value of the first exponential
    order1 = rl * il / (rd - rl) * Math.exp(-rl * t);
    // calculate the value of the second exponential
    order2 = (id - rl * il / (rd - rl)) * Math.exp(-rd * t);
    // scale it and raise the floor up to the baseline
    pred = scale * (order1 + order2) + baseline;
    // calculate the squared error and multiply by the weight (already normalized)
    error += Math.pow(pred - y, 2) * weights[index];
  }
  
  // take the square root of the error (last step of the quadratic mean)
  return Math.pow(error, 0.5)
}


/**
 * Dumb stochastic optimization algorithm. Randomly tries values around the best found, when it fails to improve it tries closer to the best and when it succeeds it tries farther. It uses the compute_weights function to pre-compute weights that compensate for non-uniform sample density along the first axis.
 *
 * @param {array[n]} x0 The initial input for the optimization (length n).
 * @param {matrix[n,2]} bounds The bounds of each variable (nx2).
 * @param {any} fun_args Fixed variables to pass to the function.  (length n).
 * @param {number} iters The number of iterations.
 * @param {boolean} maximize Whether to maximize or minimize.
 * @param {array[n]} alpha0 The initial maximum step value. (length n).
 * @param {number} gamma_better How much to increase the alpha on improvement (float).
 * @param {number} gamma_worse How much to decrease the alpha on failure (float).
 * @param {number} alpha_tol Alpha value that finishes the optimization: if the alphas reduce by more than this simulation finishes .
 * @param {function(DON'T USE)} fun The function to optimize: DO NOT USE IN THE SHEET, ONLY IN THE APP SCRIPT
 * @return {array[n+3]} result The best input found followed by the output for that input, the number of iters run, and the alpha at the end (length n+3). 
 * @customfunction
 */
function dumb_optim(x0, bounds = null, fun_args = null, iters = 100, maximize = true, alpha0 = null, gamma_better = 1.1, gamma_worse = 0.9, alpha_tol = 1e-3, fun = ldx_dx_model_error) {
  let len = x0.length;
  let alpha = (alpha0 == null) ? new Array(len).fill(1) : [...alpha0.flat()];
  let tol_tracker = 1;

  // set x0 as the best and force it to be within the bounds 
  let best_x = [...x0];
  if (bounds != null) {
    for (let i = 0; i < len; i++) {
      best_x[i] = Math.min(Math.max(best_x[i], bounds[i][0]), bounds[i][1]);
    }
  }

  // pre-compute the weights along the first axis (need to make this more modular to allow for other weighting schemes)
  full_t = fun_args.map((t) => t[0]);
  weights = compute_weights(full_t);
  // add the weights to the list of arguments to pass to the function being optimized
  fun_args = [fun_args,weights];

  // compute the error for the inital guess and start the iteration loop
  let best_value = fun(best_x.flat(), fun_args);
  for (var iter = 0; iter < iters; iter++) {
    // set the current guess to the best guess and randomly shift it (scaling with alpha)
    current_x = [...best_x];
    for (let i = 0; i < len; i++) {
      current_x[i] += (Math.random() - 0.5) * alpha[i];
    }
    // force the current guess to be within the bounds
    if (bounds != null) {
      for (let i = 0; i < len; i++) {
        current_x[i] = Math.min(Math.max(current_x[i], bounds[i][0]), bounds[i][1]);
      }
    }

    // compute the error for the current guess and check if it's the new best guess
    current_value = fun(current_x.flat(), fun_args);
    better = maximize ? (current_value > best_value) : (current_value < best_value);
    if (better) {
      // set the best guess to the current guess 
      best_value = current_value;
      best_x = [...current_x];
      // increase the alpha
      alpha = alpha.map(x => x * gamma_better);
      // update the alpha tolerance tracker
      tol_tracker *= gamma_better;
    }
    else {
      // decrease the alpha
      alpha = alpha.map(x => x * gamma_worse);
      // update the alpha tolerance tracker
      tol_tracker *= gamma_worse;
    }

    // break the loop if the alpha tolerance has been reached
    if (tol_tracker <= alpha_tol) {
      console.log("REACHED ALPHA TOL")
      break;
    }
  }

  // add the extra info to the best_x array and return it
  best_x.push(best_value);
  best_x.push(iter);
  best_x.push(tol_tracker);
  return best_x
}