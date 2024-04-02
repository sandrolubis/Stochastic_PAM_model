## A simple stochastic prototype for Propagating Annular Modes (Lubis and Hassanzadeh, 2021)

This repository contains scripts for developing a basic stochastic model for both propagating and non-propagating annular modes, drawing upon the methodologies outlined by Lubis and Hassanzadeh et al.,2021 (Eqs. 11 and 12), published in the Journal of Atmospheric Sciences (JAS).


We begin by constructing a simple stochastic system to produce synthetic time series \( z \) and \( m \) in the presence or absence of cross-EOF feedbacks. The equations of this system are the same as Eqs. (4) and (5) and (7) and (8). Following Simpson et al. (2013), we generate a synthetic time series of the random component of the eddy forcing \( \tilde{m}_{1,2} \) using a second-order autoregressive (AR(2)) noise process:

\[
\tilde{m}_1(t) = 0.6\tilde{m}_1(t - 2) - 0.3\tilde{m}_1(t - 1) + \epsilon_1(t),
\]

\[
\tilde{m}_2(t) = 0.6\tilde{m}_2(t - 2) - 0.3\tilde{m}_2(t - 1) + \epsilon_2(t),
\]

where \( t \) denotes time (in days) and \( \epsilon \) is white noise distributed uniformly between \(-1\) and \(+1\).


