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
       'postscf': 'gw',
       'selfenergy_state_range': 0,
       'xyz_file': '../h2o.xyz'
     }


print("gamma (bohr^-1)  Kohn-Sham HOMO (eV) GW HOMO (eV)")

gammas = np.arange(0., 1.01, 0.10)
egks_gamma = np.zeros(len(gammas))
egw_gamma = np.zeros(len(gammas))

for i, g in enumerate(gammas):

     ip["gamma_hybrid"]  = float(g)

     inp = Molgw_input(ip)
     run = inp.run(tmp="tmp1", keep_tmp=False)

     egks = run.homo_energy("gks")
     egw = run.homo_energy("gw")
     egks_gamma[i] = egks
     egw_gamma[i] = egw
     print(f"{g:.2f} {egks:8.3f} {egw:8.3f}")

# Spline interpolation
gks_spline = UnivariateSpline(gammas, egks_gamma, s=0)
gw_spline = UnivariateSpline(gammas, egw_gamma, s=0)
diff_spline = UnivariateSpline(gammas, egks_gamma - egw_gamma, s=0)

# Find the crossing point
gamma0 = diff_spline.roots()[0]

# Plot
gammas_fine = np.linspace(gammas[0], gammas[-1], 100)
plt.plot(gammas, egks_gamma, 'o', label='gKS', color='blue')
plt.plot(gammas_fine, gks_spline(gammas_fine), '-', color='blue')
plt.plot(gammas, egw_gamma, 'o', label='$GW$', color='green')
plt.plot(gammas_fine, gw_spline(gammas_fine), '-', color='green')
plt.plot(gamma0, gw_spline(gamma0), 's', label=f'Optimal $\\gamma$={gamma0:.3f}', color='red')
plt.xlabel("$\\gamma$ (bohr$^{-1}$)")
plt.ylabel("HOMO energy (eV)")
plt.legend()
plt.show()
