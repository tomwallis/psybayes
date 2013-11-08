\name{stan_sample}
\alias{stan_sample}
\title{Wrapper function for fitting stan models.}
\usage{
  stan_sample(file_path = NULL, data = list(),
    model_string = NULL, fit = NULL,
    n_saved_samples = 1000, iter = 1000,
    warmup = floor(iter/2), num_chains = 4,
    thin = floor((iter - warmup)/(n_saved_samples/num_chains)),
    seed = sample.int(.Machine$integer.max, 1),
    parallel = TRUE, n_cores = 4, ...)
}
\arguments{
  \item{file_path}{A character string file name or
  connection that R supports containing the text of a model
  specification in the Stan language, if model_string not
  provided.}

  \item{data}{A list object containing named data that the
  stan model expects.}

  \item{model_string}{A character vector of a stan model
  string (beginning with data block, etc.), if file not
  provided.}

  \item{fit}{A Stan fit object to take new samples from.
  This will cause file_path and model_string to be ignored}

  \item{n_saved_samples}{The number of samples to save,
  across all chains.}

  \item{iter}{The number of iterations to run the sampler
  per chain (must be at least n_saved_samples/n_chains).}

  \item{warmup}{The number of iterations that should be
  warmup.}

  \item{num_chains}{The number of independent MCMC chains
  to run.}

  \item{thin}{Stan will save every \code{thin} sample.}

  \item{seed}{The integer seed to use for seeding the
  chains.}

  \item{parallel}{Whether to run the model in parallel (one
  chain per core; uses doMC package).}

  \item{n_cores}{The number of cores to use on the machine
  you're running.}

  \item{...}{Additional arguments passed to \code{stan}
  function.}
}
\value{
  Returns an rstan fit object.
}
\description{
  A wrapper function for fitting Stan models using Tom's
  common settings. Also includes some ability to
  paralellise code. Requires the rstan package to run.
}
\examples{
data(cars)
model <- '
data{
int N;
vector[N] distance;
vector[N] speed;
}

parameters{
vector[2] beta;
real<lower=0> sigma;
}

model{
vector[N] mu;
beta ~ normal(0,100);
sigma ~ gamma(2,1e-5);
mu <- beta[1] + beta[2] * speed;
distance ~ normal(mu, sigma);
}
'

data <- list(
  N = nrow(cars) ,
  distance = cars$dist ,
  speed = cars$speed
)

# fit a linear model predicting distance from speed
fit <- stan_sample(model_string=model, data=data, iter = 250)
print(fit)

# compare to lm:
summary(lm(dist ~ speed, cars))
}
\author{
  Thomas Wallis.
}
\seealso{
  \link{stan}
}

