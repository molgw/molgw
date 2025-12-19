import numpy as np
import molgw
from molgw import Molgw_input, Molgw_output, Molgw_output_collection
from scipy.interpolate import UnivariateSpline
import matplotlib.pyplot as plt

ip = { 'scf': 'RSH',
       'basis': 'Def2-TZVPP',
       'auxil_basis': 'auto',
       'alpha_hybrid': 0.25,
       'beta_hybrid': 0.75,
       'gamma_hybrid': 1.00,
       'postscf': '',
       'selfenergy_state_range': 0,
       'xyz_file': '../h2o.xyz'
     }


print("gamma (bohr^-1)  Kohn-Sham HOMO (eV)  ΔE = E(0) - E(+) (eV)")

gammas = np.arange(0., 1.01, 0.10)
egks_gamma = np.zeros(len(gammas))
e0_gamma = np.zeros(len(gammas))
e1_gamma = np.zeros(len(gammas))

for i, g in enumerate(gammas):

     ip["gamma_hybrid"]  = float(g)

     # Neutral calculation
     ip["charge"] = 0.0
     ip["nspin"] =  1
     ip["magnetization"] = 0.0
 
     inp  = Molgw_input(ip)
     run0 = inp.run(tmp="tmp0", keep_tmp=True)
 
     egks = run0.homo_energy("gks")
     e0   = run0["scf energy"]["total"] * molgw.Ha_eV
     egks_gamma[i] = egks
     e0_gamma[i]   = e0

     # Cation calculation
     ip["charge"] = 1.0
     ip["nspin"] =  2
     ip["magnetization"] = 1.0

 
     inp  = Molgw_input(ip)
     run1 = inp.run(tmp="tmp1", keep_tmp=True)
 
     e1 = run1["scf energy"]["total"] * molgw.Ha_eV
     e1_gamma[i] = e1

     print(f"{g:.2f} {egks:8.3f} {e0-e1:8.3f}")

# Spline interpolation
deltae_gamma = e0_gamma - e1_gamma
gks_spline    = UnivariateSpline(gammas, egks_gamma, s=0)
deltae_spline = UnivariateSpline(gammas, deltae_gamma, s=0)
diff_spline = UnivariateSpline(gammas, egks_gamma - deltae_gamma, s=0)

# Find the crossing point
gamma0 = diff_spline.roots()[0]

# Plot
gammas_fine = np.linspace(gammas[0], gammas[-1], 100)
plt.plot(gammas, egks_gamma, 'o', label='gKS HOMO', color='blue')
plt.plot(gammas_fine, gks_spline(gammas_fine), '-', color='blue')
plt.plot(gammas, deltae_gamma, 'o', label='$\\Delta E$', color='green')
plt.plot(gammas_fine, deltae_spline(gammas_fine), '-', color='green')
plt.plot(gamma0, deltae_spline(gamma0), 's', label=f'Optimal $\\gamma$={gamma0:.3f}', color='red')
plt.xlabel("$\\gamma$ (bohr$^{-1}$)")
plt.ylabel("HOMO energy (eV)")
plt.legend()
plt.show()
