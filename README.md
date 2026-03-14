# Rubber Band Simulation

Monte Carlo implementation for a simulation of a 1D entropic rubber band. 

## Overview

The fundamental principle behind all Monte Carlo (MC) methods is repeated random sampling, the source and/or ''trueness'' of this randomness may vary, in order to test and find solutions to otherwise deterministic problems. MC has a vast range of applications, but in this project it was used in a relatively simple physical context where rubber bands being pulled with a force along 1 dimension were simulated. Then, by analysing how these rubber bands stretched they could be quantitatively compared with known analytical predictions and the efficiency/limitations of using MC methods could be seen. This was achieved in three parts: firstly, the MC sampling was set up by comparing an unbiased sample of rubber bands, meaning without any force, and comparing it with the known analytical results; secondly, an external force was applied via reweighing, meaning a Boltzmann weight was applied to each microstate in a large unbiased ensemble, i.e. the bias was applied after sampling; then finally the third task was to instead apply the forces directly and create the different configurations using the biased distributions. 


The simulation models a rubber band as a collection of N independent segments (links), where each segment can can be in one of two states; turned in either direction +a or -a. The behavior is random but can be influenced by an external force following statistical methods and Boltzman distribution. 


### Hooke's Law
The simulated rubber band can be greatly simplified to a chain of $N$ links along one dimension each with a length of $a$, where each link either points in the positive or the negative direction, the total length of the chain can then be defined as

$$L = a(2n - N),$$

where $n$ is the number of links pointing in the positive direction. From the assignment manual (\cite{labmanual}), it was shown that the relation between entropy and force at a fixed energy together with the relation between entropy and accessible microstates, given by the Binomial coefficient, a final expression for the relation between force $f$ and $L$ is

$$f = \frac{k_B T}{Na^2}L,$$

where $k_B$ is Boltzmann's constant, and $T$ is the temperature of the system. This expression can then be used to find Hooke's law, which states that

$$f = kL,$$

given an entropic spring constant

$$k = \frac{k_B T}{Na^2},$$

based on the parameters of the system.

### Statistical Equations
As the nature of this project is probabilistic, it is necessary to define the equations dictating the configuration and distribution of the system. Firstly, the number of microstates $\Omega$ that the rubber band can have with $n$ links in the positive direction is found using the Binomial coefficient

$$\Omega(N, n) = {N \choose n} = \frac{N!}{n!(N-n)!}.$$

Then the analytical form of the probability to find a macrostate with total length $L$ is given as

$$(L) = \frac{\Omega(N, n)}{2^N},$$

where $2^N$ is the total number of possible microstates. 

In the second and third parts, biases were applied either before or after sampling, as forces were included. The physical motivation is given by considering a single microstate with total energy $E$ and constant internal energy $E_0$. Then, by applying a constant force $f$ along the positive direction, the microstate will do work equal to $-fL$ meaning $E = E_0 - fL$. Given this $E$ and $f$, together with some temperature $T$, the Boltzmann weight, the significance, of that microstate can be described with


$$w \propto e^{\beta E} = e^{-\beta E_0}e^{\beta fL},$$

where $\beta = (k_B T)^{-1}$. In this project, $E_0$ is assumed to be zero and the constant between $w$ and $e^{\beta fL}$ is set to $1$, meaning it can be reduced to just

$$w(\text{microstate}) = e^{\beta f L}.$$

Then, the reweighed analytical probability, given total length $L$ and a force $f$, was given to be

$$P_f(L) = \frac{\Omega(N, n)e^{\beta f L}}{Z(f)}, \quad Z(f) = \sum_L\Omega(N, n)e^{\beta f L}.$$

In order to quantify how well the biased and unbiased ensembles overlap in the second part, the number of effective samples was given as

$$\mu_\text{eff} = \frac{1}{\sum_i \bar w_i^2},$$

where

$$\bar w_i = \frac{w_i}{\sum_j w_j}.$$

$\mu_\text{eff}$ can test the efficiency of the reweighing with $M$ samples, and if $\mu_\text{eff} \approx M$, then there is good overlap, but if $\mu_\text{eff} \ll M$, there is bad overlap and resampling might be required. 

## Features

- **Monte Carlo implementation**: Implements a simple MC method to simulate a rubber band with statistical physics background
- **Calculate probabilities for chain lenghts**: Gives the probability for a given chain length with given number of links for various external forces
- **Investigate efficiency of reweighting**: Examines the efficiency of reweighting for increasing applied forces


## Requirements
```
numpy
matplotlib
scipy
math
os
```

Install dependencies with:
```bash
pip install -r requirements.txt
```

## Usage
We set a specific seed for reproducability:
```python
rng = np.random.default_rng(42)
```
Then we define the classes and functions that help with calculating the weighting efficiency and the theoretical probability.

### Unbiased rubber bands

First, we considered an ensemble of unbiased microstates: rubber bands which are not pulling with any force and are just random configurations. We set up a `Link` and `RubberBand` class to simulate the rubber bands and links between them. We calculated the length of a rubber band. We did a Monte Carlo simulation of $M = 10^5$ unbiased rubber bands with $N = 100$ links and $a = 1$. The exact probability distribution along with the Monte Carlo samples and the ratio of the Monte Carlo and the exact solution can be seen in the figure below. 

![Length distribution of an unbiased rubber band. Each rubber band has $N = 100$ links with $a = 1$, and $M = 10^5$ Monte Carlo samples were made. The Monte Carlo histogram is overlaid with the exact solution. The ratio plot shows good agreement with the analytical solution, along with the $\chi^2/$ndf value of $0.745 \sim 1$.](figures/I_ratio.pdf)

We can see an agreement between the Monte Carlo sampling and the exact analytical model. The ratio plot in the lower panel of the figure demonstrates that the model is in good agreement with the MC data (although it is seen to be diverging beyond the region of statistical convergence). Furthermore, the $\chi^2/$ndf shows a value of $0.745 \sim 1$, which also shows good agreement. 

### Reweighted pulling with a force

Now, we will introduce an external force $f$ pulling the rubber band in the positive direction. We will do this with a reweighting technique, such that we do not need to resample the rubber bands, but rather we calculate a new weight with which they enter. We will add a Boltzmann weight to each `RubberBand`. We will work in units where $k_B = 1$ and $T = 1$, such that $\beta = 1$. We compare the weighted analytical solution for different forces, $f = \{0.01, 0.05, 0.1, 0.5, 1.0\}$, which can be seen in the figures below. For low external forces, the analytical solution is quite close to the MC sample, while for higher forces (from $\sim 0.5$) the analytical distribution diverging from the MC samples. \par 

![Reweighted length distributions $P_f(L)$ for various external forces; $f \in \{0.01, 0.05, 0.1, 0.5, 1.0\}$. For low external forces, the analytical solution is quite close to the MC sample, while for higher forces (from $\sim 0.5$) the analytical distribution diverging from the MC samples. ](figures/II_f=0.01.pdf)

![$f = 0.05$ ](figures/II_f=0.05.pdf)

![$f = 0.1$ ](figures/II_f=0.1.pdf)

![$f = 0.5$ ](figures/II_f=0.5.pdf)

![$f = 1.0$ ](figures/II_f=1.0.pdf)


We can see the reweighting efficiency in the figure below. A rapid decay can be seen for higher forces, which indicates that the reweighting becomes statistically unreliable. This is likely because only a tiny fraction of the original random configurations have a sufficiently large length to contribute significant weight under higher forces. This causes the statistical error to grow, and therefore the ratio plot to diverge. 

![Reweighting efficiency $(\mu_{\rm eff}/M)$ as s function of the applied force $f$. A rapid decay can be seen for higher forces, which indicates that the reweighting becomes statistically unreliable.](figures/II_efficiency.pdf)


### Direct simulation of pulling with a force

We continue our simulation of the rubber band with a constant applied force $f$, now drawing Monte Carlo samples using a biased distribution to obtain more reliable results for high $f$. We modify our `RubberBand` class to take both `kbT` and `force` inputs --- the former is set to $1.0$ for our experiments, while the latter are sampled from the distribution $f\in[0.0,2.5]$. When `force` is given, we draw string directions from `numpy.rng.choice` with a biased `p`-parameter:

$$p=\{p_+, \; p_- \} = \{p_+,\;1-p_+\}$$

$$p_+=\frac{1}{2}\left(1+\tanh(\beta fa)\right)$$

For comparison against expectation, we initialize $5000$ rubber bands at each value of `force`. Forces are drawn from the array $(f\in[0.01, 0.05],\;\texttt{step}=0.001)\cup(f\in[0.5, 1.0],\;\texttt{step}=0.01)\cup(f\in[1.0, 2.5],\;\texttt{step}=0.1)$. These bins mean we have the highest precision in the linear regime at low $f$. The mean extension of each rubber band is then simply the arithmetic mean of the extension from all the generated rubber bands (in our code, this is the mean of the ensemble of \texttt{RubberBand.length()} generated for each force). We then also calculate the standard deviation of these values, which we use as our error during comparison versus analytical models.

We get an analytical function $\left<L\right>(f)=Na\tanh(\beta fa)$, which has a roughly linear regime for low $f$ and then transitions to a constant regime for high $f$. In the figure below, this is plotted as the red line, while the results from our Monte-Carlo sampling are the blue dots. The standard deviation of the Monte-Carlo samples are visualized as a transparent area around the main function.

![Comparison of Monte-Carlo and analytical results. With the full expanded analytical mean $\left<L\right>$, we find a very good agreement for all values of $f$. The standard deviation is initially quite large at low values of force, and gets narrower over time.](figures/III_mean_L_force_0.pdf)

Defining the ''linear regime'' as the region between $0.01$ and $0.05$, we can perform a linear fit to our Monte-Carlo sample pulled only from forces in this range (conveniently, this is the first element of our force array), using `SciPy`'s module `curve_fit`. Note that in the figure below, we take the domain of the plot as the range of the first two elements of the force array, $(0.0,\;0.5]$ --- but only the first component $(0.0,\;0.05]$ is fitted.

![Linear-regime fit for $\left<L\right>$. The analytical value for the constant $k$ is extremely close to the fitted value --- the absolute difference between $k=\frac{k_BT}{Na^2}$ and our calculated fit for $k_\textup{eff}$ is $\approx2.47\times10^{-5}$. If we instead fit the entire range of the plot, the absolute difference becomes $\approx2.36\times10^{-3}$. (These values are sensitive to seeding.)](figures/III_mean_L_force_1.pdf)

For large $f$, the linear fit breaks down, visualized in the figure below. The distribution becomes much tighter --- standard deviations drop by an order of magnitude from the linear regime to the large-force regime (from $\sigma_\text{linear}\approx10.18$ to $\sigma_\textup{large}\approx1.81$). This is expected from analytical results as well.

![Linear-regime fit for $\left<L\right>$ with the plot extended to higher forces. The disagreement between the linear regime fit and the final Monte Carlo / expanded analytical derivations is roughly $\Delta\left<L\right>\approx141.62$, while the initial value for this quantity is $\Delta\left<L\right>\approx0.018$, a difference of five orders of magnitude. The main body of divergence starts after $f\sim0.5$, which has a disagreement $\Delta\left<L\right>\approx3.80$. (These values are also sensitive to seeding.)](figures/III_mean_L_force_2.pdf)


## Classes

### `Link`

Represents an individual link segment with two states.

**Parameters:**
- `direction` (int, +-1): direction of the link

### `RubberBand`

Simulates a chain of N link segments.

**Parameters:**
- `N` (int): Number of links in the chain. Default: 100
- `a` (float, default 1.0): Link object defining segment properties
- `kBT` (float, default None): Thermal energy (Boltzmann constant × Temperature). Default=1
- `force` (float, default None): Force
- `rng` (rng, default None): Random number generator, if set to None `np.random.default_rng(seed=42)` will be used.
If the bias and kbT are given, `RubberBand` will create a biased sample. If not, the `RubberBand` will be unbiased. 

**Methods:**
- `length()`: Calculate total length of the chain
- `weigt(force, T, kb=1.0)`: Compute Boltzmann weight for given force


## Project Structure
```
.
├── entropic_spring.py # Main simulation code
└── README.md         # This file
└── requirements.txt # required modules for the code
└── Figures     # All produced Figures
```















\section{Discussion and Conclusion}
Our $\chi^2/\text{ndf}$ value for the unbiased rubber band showed great agreement with the theoretical distribution, meaning that our sampling was indeed both efficient enough and convergent enough to accurately model the system. For smaller values of $M$, we found that the $\chi^2/\text{ndf}$ value was lower, as expected (around $0.745$ for $M=1\times10^5$ samples). (Larger values of $M$ were not possible due to computational constraints.) For the case of applied constant forces and reweighted rubber bands, we found a decent agreement with the analytical distribution for low forces, but one that breaks down in the case of large forces, as expected due to the unbiased and biased ensembles overlapping less. This leads directly into our true simulation of the biased distribution, which was able to correct this artifact, and more effectively model the system. Our results accurately recreate both the full expected analytical distribution ($\tanh$) and the Taylor-expanded linear regime, with an extremely accurate derived value for the effective Hooke's constant. We once again find the very considerable divergence (several orders of magnitude with $\beta=\frac{1}{k_BT}=1$) for large forces (larger than $f=0.5$), which means the Hooke's law approximation should break down in similar regimes when either $\beta$ or $f$ are large. This has implications especially for low-energy physics, where $T$ is small and models involving Hooke's law as an approximation may begin to break down.
