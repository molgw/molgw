import numpy as np
from . import __version__, Ha_eV


def evaluate(input_file="dipole_time.dat", output_file="spectrum.dat",
             npadding=10, damping_time=None, ntruncate=None ):
    """
       Calculate the optical spectrum in frequency from a dipole series in time
    """
    
    # Read RT data
    dipole = np.genfromtxt(input_file)
    print(f"Opening file: {input_file}")
    nt_read = len(dipole[:,0])
    print(f"Number of time steps found: {nt_read}")

    # Tuning ntruncate can cut the end of the trajectory 
    if ntruncate is None:
        nt = nt_read
    else:
        if ntruncate < nt_read:
            print(f'Truncate the time series after snapshot {ntruncate}')
            nt = min(nt_read, ntruncate)
        else:
            nt = nt_read
    
    T    = dipole[:nt,0]
    dt   = T[1] - T[0]
    tsim = T[nt-1] - T[0]
    
    # Tuning damping_time or take it as 0.1 * tsim
    if damping_time is None:
        tau = tsim * 0.1
    else:
        tau = damping_time

    eta = 1.0 / tau * Ha_eV
    print(f"Damping time: {tau:.2f} (a.u.)")
    print(f"that corresponds to a spectrum broadening: {eta:.4f} (eV)")
    
    
    print(f"Zero padding with factor: {npadding}")
    nt_with_padding = nt * npadding
    
    Ex   = dipole[:nt,1]
    Ey   = dipole[:nt,2]
    Ez   = dipole[:nt,3]
    Dx   = dipole[:nt,4] - dipole[0,4]
    Dy   = dipole[:nt,5] - dipole[0,5]
    Dz   = dipole[:nt,6] - dipole[0,6]
    
    #  Calculate the discrete energies in Ha and eV
    Omega = np.zeros((nt_with_padding))
    for i in range(npadding*nt):
        Omega[i] = 2.0 * np.pi * i / (tsim * npadding)
    Omega_eV = Omega * Ha_eV
    
    ########################
    # Damping
    Dx *= np.exp( -T / tau )
    Dy *= np.exp( -T / tau )
    Dz *= np.exp( -T / tau )
    
    ########################
    # Zero padding (Fourier interpolation)
    Ex = np.concatenate( (Ex, np.zeros(((npadding-1)*nt))) )
    Ey = np.concatenate( (Ey, np.zeros(((npadding-1)*nt))) )
    Ez = np.concatenate( (Ez, np.zeros(((npadding-1)*nt))) )
    Dx = np.concatenate( (Dx, np.zeros(((npadding-1)*nt))) )
    Dy = np.concatenate( (Dy, np.zeros(((npadding-1)*nt))) )
    Dz = np.concatenate( (Dz, np.zeros(((npadding-1)*nt))) )
    
    # Fourier transforms
    Dx_FFT = np.fft.fft(Dx, norm='forward')
    Dy_FFT = np.fft.fft(Dy, norm='forward')
    Dz_FFT = np.fft.fft(Dz, norm='forward')
    
    Ex_FFT = np.fft.fft(Ex, norm='forward')
    Ey_FFT = np.fft.fft(Ey, norm='forward')
    Ez_FFT = np.fft.fft(Ez, norm='forward')
    
    # Spectrum: sum of the imaginary parts in x, y, z directions
    spectrum = np.imag(Dx_FFT / Ex_FFT) + np.imag(Dy_FFT / Ey_FFT) + np.imag(Dz_FFT / Ez_FFT)
    # FIXME Should we do that?
    spectrum *= Omega
    
    nomega_plot = int(nt_with_padding / 2)
    
    header = "Frequency (eV)   Spectrum (a.u.)"
    data = np.column_stack((Omega_eV[:nomega_plot], spectrum[:nomega_plot]))
    np.savetxt(output_file, data, delimiter='    ', header=header, fmt='% .5e')

    return Omega_eV[:nomega_plot], spectrum[:nomega_plot]

