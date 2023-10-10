import pycwt
import numpy as np

def pad(x, method='reflection'):
    """Pad to bring the total length N up to the next-higher 
    power of two.

    From David Albanese mlpy

    :Parameters:
       x : 1d array_like object
          data
       method : string ('reflection', 'periodic', 'zeros')
          method
          
    :Returns:
       xp, orig : 1d numpy array, 1d numpy array bool
          padded version of `x` and a boolean array with
          value True where xp contains the original data
    """

    x_arr = np.asarray(x)

    if not method in ['reflection', 'periodic', 'zeros']:
        raise ValueError('method %s not available' % method)
    
    diff = next_p2(x_arr.shape[0]) - x_arr.shape[0]
    ldiff = int(diff / 2)
    rdiff = diff - ldiff

    if method == 'reflection':
        left_x = x_arr[:ldiff][::-1]
        right_x = x_arr[-rdiff:][::-1]         
    elif method == 'periodic':
        left_x = x_arr[:ldiff]
        right_x = x_arr[-rdiff:]
    elif method == 'zeros':
        left_x = np.zeros(ldiff, dtype=x_arr.dtype)
        right_x = np.zeros(rdiff, dtype=x_arr.dtype)
        
    xp = np.concatenate((left_x, x_arr, right_x))
    orig = np.ones(x_arr.shape[0] + diff, dtype=np.bool)
    orig[:ldiff] = False
    orig[-rdiff:] = False
   
    return xp, orig


def autoscales(N, dt, dj, wf, p):
    """Compute scales as fractional power of two.
    :Parameters:
       N : integer
          number of data samples
       dt : float
          time step
       dj : float
          scale resolution (smaller values of dj give finer resolution)
       wf : string
          wavelet function ('morlet', 'paul', 'dog')
       p : float
          omega0 ('morlet') or order ('paul', 'dog')

    :Returns:
       scales : 1d numpy array
          scales
    """
    if wf == 'dog':
        s0 = (dt * np.sqrt(p + 0.5)) / np.pi
    elif wf == 'paul':
        s0 = (dt * ((2 * p) + 1)) / (2 * np.pi)
    elif wf == 'morlet':
        s0 = (dt * (p + np.sqrt(2 + p**2))) / (2 * np.pi)
    else:
        raise ValueError('wavelet function not available')
    J = int(np.floor(int(dj**-1) * np.log2((N * dt) / s0)))
    s = np.empty(J + 1)
    for i in range(s.shape[0]):
        s[i] = s0 * 2**(i * dj)
    return s

def next_p2(n):
    """Returns the smallest integer, greater than n
    (n positive and >= 1) which can be obtained as power of 2.
    """
    
    if n < 1:
        raise ValueError("n must be >= 1")
    v = 2
    while v <= n:
        v = v * 2
    return v


def recon(wf, p, dj, dt):
    if wf == 'morlet':
        if p != 6:
            print('ERROR P NOT EQUAL TO 6, ICWT WILL LOOK INCORRECT. CONTINUING ANYWAY...')
        rf = (dj*(dt**0.5))/(0.776*(np.pi**-0.25))  # morlet: equation 11 in Torrance and Compo
    if wf == 'dog':
        if p == 2:
            rf = (dj*(dt**0.5))/(3.541*0.867)  # dog: equation 11 in Torrance and Compo
        if p == 6:
            rf = (dj*(dt**0.5))/(1.966*0.884)
        if p != 2 and p !=6:
            print('ERROR P NOT EQUAL TO 2 or 6, ICWT WILL LOOK INCORRECT. CONTINUING ANYWAY...')
            rf = (dj*(dt**0.5))/(3.541*0.867)
    return rf


def fourier_from_scales(scales, wf, p):
    """Compute the equivalent Fourier period
    from scales.
    
    :Parameters:
       scales : list or 1d numpy array
          scales
       wf : string ('morlet', 'paul', 'dog')
          wavelet function
       p : float
          wavelet function parameter ('omega0' for morlet, 'm' for paul
          and dog)
    
    :Returns:
       fourier wavelengths
    """

    scales_arr = np.asarray(scales)

    if wf == 'dog':
        return  (2 * np.pi * scales_arr) / np.sqrt(p + 0.5)
    elif wf == 'paul':
        return  (4 * np.pi * scales_arr) / float((2 * p) + 1)
    elif wf == 'morlet':
        return  (4 * np.pi * scales_arr) / (p + np.sqrt(2 + p**2))
    else:
        raise ValueError('wavelet function not available')

def icwt_fixed(W, sj, dt, dj=1/12, wavelet='morlet'):
    """Inverse continuous wavelet transform.
    Parameters
    ----------
    W : numpy.ndarray
        Wavelet transform, the result of the `cwt` function.
    sj : numpy.ndarray
        Vector of scale indices as returned by the `cwt` function.
    dt : float
        Sample spacing.
    dj : float, optional
        Spacing between discrete scales as used in the `cwt`
        function. Default value is 0.25.
    wavelet : instance of Wavelet class, or string
        Mother wavelet class. Default is Morlet
    Returns
    -------
    iW : numpy.ndarray
        Inverse wavelet transform.
    Example
    -------
    >> mother = wavelet.Morlet()
    >> wave, scales, freqs, coi, fft, fftfreqs = wavelet.cwt(var,
           0.25, 0.25, 0.5, 28, mother)
    >> iwave = wavelet.icwt(wave, scales, 0.25, 0.25, mother)
    """
    wavelet = pycwt.wavelet._check_parameter_wavelet(wavelet)

    a, b = W.shape
    c = sj.size
    if a == c:
        sj = (np.ones([b, 1]) * sj).transpose()
    elif b == c:
        sj = np.ones([a, 1]) * sj
    else:
        raise Warning('Input array dimensions do not match.')

    # As of Torrence and Compo (1998), eq. (11)
    iW = (dj * np.sqrt(dt) / (wavelet.cdelta * wavelet.psi(0)) * (np.real(W) / np.sqrt(sj)).sum(axis=0))

    return iW







