psybayes
========

An R package for modelling psychophysical data using MCMC.

psybayes is a package to facilitate the fitting of data from psychophysical experiments in a Bayesian framework.
Models are fit using MCMC sampling conducted through Stan (and the rstan package).
These will need to be properly installed on your system first (see http://www.mc-stan.org for setup).
Many of these functions are wrappers for rstan calls, or plotting functions that are commonly used.

Todo
========

  * Implement WAIC in binomial_model_comparison_metrics.
  * Bring mcmc helper functions into package, with documentation.
  * Bayesian correlation computation functions.


License
========

This is free, open source software released under the GPL-3 License. See LICENSE.txt for details.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE
LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

Use this software at your own risk.
