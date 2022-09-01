def dRdu(R):
    """
    Parameters:
    -----------
    R: array
        Values of resistances.
    
    Returns:
    --------
    derivative: array
        The numerical derivative of R, logarithmically transformed.
    """
    log = np.log(R)
    derivative = np.gradient(log)
    return derivative
    
def Mott(f, T):
    """
    Mott formula to calculate the Seebeck coefficient
    
    Parameters:
    -----------
    f: array
        The numerical derivative of log(R).
    T: array
        The temperature range for the Seebeck coefficient.
        
    Returns:
    --------
    mott: array
        An array with Seebeck coefficients for different temperature values.
    """
    mott = pi**2 * k_b**2 * T / (3*e) * f
    return mott

def EfBLG(Vg, Vg0):
    """
    Models Fermi energy for 2D materials for a given series of voltages
    
    Parameters:
    -----------
    Vg: array
        Applied top gate voltages. 
    Vg0: float
        Symmetric voltage Vg0.
    
    Returns:
    --------
    Ef: array
        Fermi energy for the 2D material (Bilayer Graphene).
    R: array
        Resistance for the 2D material, given a series of top gate voltages.
    """
    gamma1 = 0.4
    Vf0 = 1e6
    er = 3.6
    tbg = 20e-9
    n = epsilon_0*er*(Vg-Vg0)/(e*tbg)
    alpha = 1/(pi*hbar_ev**2*Vf0**2)
    n0 = epsilon_0*er*(0.5-Vg0)/(e*tbg)
    n_eff = np.sqrt(n**2 + n0**2)    
    Ef = np.sign(n)*0.5*(-gamma1+np.sqrt(gamma1**2+4*np.abs(n_eff)/alpha)) # conductance band 
    R = 1/(n_eff*e*Vf0)
    return Ef, R

# for 1D
def r_QPC(Vg):
    """
    Calculates the conductance of QPC.
    
    Parameters:
    -----------
    Vg: array
        Applied top gate voltages.
    Returns:
    --------
    R: array
        Resistance for the 1D material (Quantum point contact), given a series of top gate voltages.
    """
    Vg0=0
    m_star = gamma1/(2*Vf0**2) # band structure, all parabollic
    n = epsilon_0*er*(Vg-Vg0)/(e*tbg)
    alpha = 1/(pi*hbar_ev**2*Vf0**2)
    Ef = np.sign(n)*0.5*(-gamma1+np.sqrt(gamma1**2+4*np.abs(n)/alpha))
    N = np.round(np.sqrt(2*np.abs(Ef)*m_star)*W/(pi*hbar_ev))
    sigma = 4*e**2*N/h
    resistance = 1/sigma
    return resistance

# for 0D
def QD(x, x0, gamma, scaling):
    """
    Define a series of Lorentzian peaks to simulate the 0D behavior.
    
    Parameters:
    ----------
    x : array
        range of points on the x-axis, independent variable
    x0 : list
        the centers of the Lorentzian peak
    scaling : float
        Some scaling factor
    """
    return scaling/pi * (0.5*gamma / ((x-x0)**2 + (0.5*gamma)**2))
