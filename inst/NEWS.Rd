Changes in version 1.1.3

* Updates in INLA means now the Gompertz model is also available for
survival modelling. fit.models(), make.surv() and print() have now
been updated so that the Gompertz model can be run under 'survHE'.

* Related to this, to improve computational stability, *all* the 
INLA models are now run by 'survHE' using the following trick: 
first the times are rescaled (on the fly) in the interval [0-1] 
(by simply recomputing 'time=time/max(time'). The resulting 
models are *not* directly comparable to other inferential engines
(because they are fitted to different data), but 'survHE' 
automatically rescales the estimates and model fitting statistics
(eg *ICs) so that the 'plot' and 'print' methods give the 
correct answers. make.surv() has also been updated to reflect
this.

* There are also changes to make.transition.probs(), which 
has been updated and streamlined to compute transition 
probabilities off the survival curves fitted in a 'survHE' 
object. The computation is quicker and now based on the 
more robust relationship between the cumulative hazard function
and the transition probabilities.

* A new function make_data_multi_state() to create a dataset
in the format required to analyse data in a multi-state 
framework.

* A new function three_state_mm() added to fit a standard
3-states Markov model, based on survival curves that are then
mapped onto transition probabilities (this function is under
testing, though).
