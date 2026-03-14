import numpy as np 
import math
import matplotlib.pyplot as plt 
import os 
from scipy.optimize import curve_fit

# to ensure there is a ./plots/ directory to save the plots
plots_dir = './figures'
os.makedirs(plots_dir, exist_ok=True)

# seed for reproducability
rng = np.random.default_rng(42)

class Link():
    def __init__(self, direction):
        # direction is either +1 or -1
        self.direction = direction
    

class RubberBand():
    def __init__(self, N, a=1, kbT=None, force=None, rng=None):
        # setting up rng 
        if rng is None:
            self.rng = np.random.default_rng(seed=42)
        else:
            self.rng = rng 

        # setting up other parameters
        self.N = N
        self.a = a 

        # getting random directions; use bias if kbT and force are given
        if not(kbT and force):
            dirs = self.rng.choice([-1, 1], size=N)
        else:
            p_pos = 0.5*(1+np.tanh((1/kbT)*force*a))
            p_neg = 1 - p_pos
            dirs = self.rng.choice([-1, 1], size=N, p=[p_neg, p_pos])

        # setting up the links
        self.links = [Link(d) for d in dirs]

    
    def length(self):
        # calculating lenght of the band
        # L = a(2n-N)
        n = sum(1 for link in self.links if link.direction == 1)
        self.len = self.a * (2 * n - self.N)
        return self.len 


    def weight(self, force, T, kb=1.0):
        # calculating the boltzmann weight
        # w = exp(beta * f * L)
        beta = 1.0/(kb * T) 
        if self.len is not None: 
            return np.exp(beta * force * self.len)
        else: 
            return np.exp(beta * force * self.length())


def analytical(N, lengths_mc, a=1):
    lengths = np.arange(-N*a, N*a + 1, 2*a)
    # getting analytical model
    n = ((lengths/a + N) / 2).astype(int)
    probs = [math.comb(N, nv) / 2**N for nv in n]
    return lengths, probs


def weighted_analytical(N, f, T, a=1, kb=1.0):
    # P(L) = Omega(N, n) * exp(beta*f*L) / Z(f)
    beta = 1 / (kb * T)
    lengths = np.arange(-N*a, N*a + 1, 2*a)
    # getting analytical model
    n = ((lengths/a + N) / 2).astype(int)

    # calculate unnormalized weights
    p_unnorm = np.array([math.comb(N, nv) * np.exp(beta * f * lv) for nv, lv in zip(n, lengths)])
    
    # normalize 
    probs = p_unnorm / np.sum(p_unnorm)
    return lengths, probs

def mean_L_analytical(N, f, kbT, a):
    # mean length analytical
    beta = 1 / kbT
    return(beta*N*a*a*f)

def mean_L_exp(N, f, kbT, a):
    # mean length exponent
    beta = 1 / kbT
    return(N*a*np.tanh(beta*f*a))


if __name__ == '__main__':
    # model parameters 
    a = 1  # link length 
    N = 100  # number of links
    M = int(1e5)  # number of samples
    rng = np.random.default_rng(seed=42)

    # MC a large number of rubber bands
    rubber_bands = [RubberBand(N=N, a=a, rng=rng) for _ in range(M)]
    lengths = [rb.length() for rb in rubber_bands]
    lengths_an, p_an = analytical(N, lengths_mc =lengths, a=a)

    # plotting
    fig, ax = plt.subplots(2, 1, figsize=(6,6),
                           sharex=True, height_ratios=[2, 1])
    
    # setting probablitly for histogram
    bin_width = lengths_an[1] - lengths_an[0]    
    edges = np.concatenate(([lengths_an[0] - bin_width/2], 
                           (lengths_an[:-1] + lengths_an[1:])/2,
                           [lengths_an[-1] + bin_width/2]))
    counts, _ = np.histogram(lengths, bins=edges, density=True)
    dcounts, _ = np.histogram(lengths, bins=edges)
    pmc = counts/counts.sum()
    dpmc = dcounts/dcounts.sum()
    
    # Make the ratio mc/th for plotting
    ratio = np.array([mc / th for mc, th in zip(dpmc,p_an)])
    mask = (ratio > 0)
    # Make a chi2 test
    chi2 = 0
    ndf = 0
    for mc, th in zip(dcounts, p_an):
        pth = len(lengths) * th # calculate chi2 on counts. rescale th value
        if pth > 5: # guard against (almost) empty bins
            chi2 += (mc - pth)**2 / pth
            ndf += 1
    ndf -= 1 # here you would subtract 1 per parameter, one for the normalization

    ax[0].step(lengths_an, pmc, where='mid', label="MC")
    ax[0].plot(lengths_an, p_an, label='exact') 
    #ax[0].set_xlabel("length (m)")   
    ax[0].set_ylabel("probability")
    ax[0].legend()
   
    # Put the chi2 on the plot
    s = f"{chi2/ndf:.3f}"
    ax[0].text(0.2,0.92, r"$\chi^2/$ndf = "+s, horizontalalignment='center',
                verticalalignment='center', transform=ax[0].transAxes)

    ax[1].axhline(1.0, color='gray', linestyle='--')
    ax[1].plot(lengths_an[mask], ratio[mask])
    ax[1].set_xlabel("length/a")
    ax[1].set_ylabel("ratio MC/exact")

    plt.suptitle(f"unbiased rubber band, N={N}, a={a}, M={M}")
    plt.savefig('figures/I_ratio.png', bbox_inches='tight', dpi=200)
    #plt.savefig('figures/I_ratio.pdf', bbox_inches='tight', dpi=200)
    #plt.show()
    plt.close()

    T = 1
    forces = [0.01, 0.05, 0.1, 0.5, 1.0]
    lengths_unb = np.array([rb.length() for rb in rubber_bands])

    mu_effs = []
    for f in forces: 
        # calculating weights
        weights = np.array([rb.weight(force = f, T=T) for rb in rubber_bands])

        # effective sample size
        weights_norm = weights / np.sum(weights)
        mu_effs.append((1.0/np.sum(weights_norm**2)) / M)

        # getting analytical result
        lengths_an, p_an = weighted_analytical(N, f, a=a, T=T)

        # plotting comparison
        fig, ax = plt.subplots(2, 1, figsize=(7,7), sharex=True, height_ratios=[2, 1])
        counts, bins, _ = ax[0].hist(lengths_unb, bins=len(lengths_an), weights=weights,
                                     density=True, label=f"MC reweighted")
        ax[0].plot(lengths_an, p_an, '-', label="weighted analytical")
        ax[0].set_ylabel('probability')
        ax[0].legend()

        p_mc, _ = np.histogram(lengths_unb, bins=bins, weights=weights, density=True)
        mask = p_an > 1e-5
        ratio = p_mc[mask] / p_an[mask]
        ax[1].plot(lengths_an[mask], ratio)
        ax[1].axhline(1.0, color='gray', linestyle='--')
        ax[1].plot(lengths_an[mask], ratio)
        ax[1].set_xlabel("length/a")
        ax[1].set_ylabel("ratio MC/exact")

        plt.suptitle(f"force f = {f}")
        plt.savefig(f'figures/II_f={f}.png', bbox_inches='tight', dpi=200)
        #plt.savefig(f'figures/II_f={f}.pdf', bbox_inches='tight', dpi=200)
        #plt.show()
        plt.close()

    # effective samples plot
    plt.figure(figsize=(6, 4))
    plt.plot(forces, mu_effs, 'o-')
    plt.xlabel('force (f)')
    plt.ylabel(r"$\mu_{\text{eff}} / M$")
    plt.title("reweighting efficiency vs force")
    plt.savefig('figures/II_efficiency.png', bbox_inches='tight', dpi=200)
    #plt.savefig('figures/II_efficiency.pdf', bbox_inches='tight', dpi=200)
    #plt.show()
    plt.close()

    # function for fitting
    def f(x, a, b):
        return(a*x+b)
    
    M = int(5e3)
    kbT = 1.0

    # create array of forces
    farr_0 = np.arange(0.01, 0.05, 0.001)
    forces_0 = farr_0.tolist()
    farr_1 = np.arange(0.05, 1.0, 0.01)
    forces_1 = farr_1.tolist()
    farr_2 = np.arange(1.0, 2.5, 0.1)
    forces_2 = farr_2.tolist()
    forces = forces_0 + forces_1 + forces_2

    # initialization
    mean_Ls = []
    sd_Ls = []
    anm_Ls = []
    ane_Ls = []

    for force in forces:
        # generating new biased rubber band 
        rubber_bands = [RubberBand(N=N, a=a, kbT=kbT, force=force, rng=rng) for _ in range(M)]
        lengths_b = np.array([rb.length() for rb in rubber_bands])
        
        # calculate mean and STD
        mean_L = 1/M*np.sum(lengths_b)
        sd_L = np.sqrt(1/M*np.sum((lengths_b-mean_L)**2))
        # theoretical values
        anm_L = mean_L_analytical(N=N, f=force, kbT=kbT, a=a)
        ane_L = mean_L_exp(N=N, f=force, kbT=kbT, a=a)

        mean_Ls.append(mean_L)
        sd_Ls.append(sd_L)
        anm_Ls.append(anm_L)
        ane_Ls.append(ane_L)

    # effective spring constant
    ann = kbT/(N*a*a)

    # fit mc results
    fitted_force = forces_0
    popt, pcov = curve_fit(f, forces[:len(fitted_force)], mean_Ls[:len(fitted_force)])
    k_eff = 1/popt[0]
    x = np.linspace(0.01, 2.5, 100)

    # printing results
    print(np.abs(ann-k_eff))
    print(f'Initial STD: {sd_Ls[0]:.3f},\nFinal STD: {sd_Ls[-1]:.3f}')
    print(f'Initial Agreement: {np.abs(anm_Ls[0]-mean_Ls[0]):.3f},\nFinal Agreement: {np.abs(anm_Ls[-1]-mean_Ls[-1]):.3f}')
    print(f'0.5 Agreement: {np.abs(anm_Ls[forces.index(0.5000000000000001)]-mean_Ls[forces.index(0.5000000000000001)]):.3f}')

    # plotting
    plt.figure(figsize=(6, 4))
    plt.plot(forces, mean_Ls, 'o', color='tab:blue', alpha=0.25, markeredgecolor='k')
    plt.plot(forces, mean_Ls, '-', color='tab:blue', label='MC Mean L')
    plt.fill_between(forces, np.array(mean_Ls)-np.array(sd_Ls), np.array(mean_Ls)+np.array(sd_Ls), color='tab:blue', alpha=0.2)
    plt.plot(forces, anm_Ls, '-', color='tab:red', label='Approximate Analytical Mean L\n'+r'$\frac{k_BT}{Na^2}=$'+f'{ann:.3f}')
    plt.plot(x, f(x, popt[0], popt[1]), '--', color='tab:green', label='MC Linear Fit\n'+r'$k_\text{eff}=$'+f'{k_eff:.3f}')
    plt.xlabel('force (f)')
    plt.ylabel(r"mean length ($\left<L\right>$)")
    plt.title("Force vs. Mean Length, Linear Regime Divergence")
    plt.legend()
    plt.savefig('figures/III_mean_L_force_2.png', bbox_inches='tight', dpi=200)
    #plt.savefig('figures/III_mean_L_force_2.pdf', bbox_inches='tight', dpi=200)
    #plt.show()
    plt.close()

    plt.figure(figsize=(6, 4))
    plt.plot(forces, mean_Ls, 'o', color='tab:blue', alpha=0.25, markeredgecolor='k')
    plt.plot(forces, mean_Ls, '-', color='tab:blue', label='MC Mean L')
    plt.fill_between(forces, np.array(mean_Ls)-np.array(sd_Ls), np.array(mean_Ls)+np.array(sd_Ls), color='tab:blue', alpha=0.2)
    plt.plot(forces, ane_Ls, '-', color='tab:red', label='Expanded Analytical Mean L\n'+r'$Na\text{tanh}(\beta fa)$')
    plt.xlabel('force (f)')
    plt.ylabel(r"mean length ($\left<L\right>$)")
    plt.title("Force vs. Mean Length, Comparison with Analytical")
    plt.legend()
    plt.savefig('figures/III_mean_L_force_0.png', bbox_inches='tight', dpi=200)
    #plt.savefig('figures/III_mean_L_force_0.pdf', bbox_inches='tight', dpi=200)
    #plt.show()
    plt.close()

    plt.figure(figsize=(6, 4))
    plt.plot(forces, mean_Ls, 'o', color='tab:blue', alpha=0.25, markeredgecolor='k')
    plt.plot(forces, mean_Ls, '-', color='tab:blue', label='MC Mean L')
    plt.fill_between(forces, np.array(mean_Ls)-np.array(sd_Ls), np.array(mean_Ls)+np.array(sd_Ls), color='tab:blue', alpha=0.2)
    plt.plot(x, f(x, popt[0], popt[1]), '--', color='tab:green', label='MC Linear Fit\n'+r'$k_\text{eff}=$'+f'{k_eff:.3f}')
    plt.xlabel('force (f)')
    plt.ylabel(r"mean length ($\left<L\right>$)")
    plt.title("Force vs. Mean Length, Linear Regime Fit")
    plt.xlim(0.0,0.5)
    plt.ylim(0.0, 50.0)
    plt.legend()
    plt.savefig('figures/III_mean_L_force_1.png', bbox_inches='tight', dpi=200)
    #plt.savefig('figures/III_mean_L_force_1.pdf', bbox_inches='tight', dpi=200)
    #plt.show()
    plt.close()