function shakuradelay(lambdas::Vector{T} where T<:Real; M=M, R=R)

    map(l->shakuradelay(l; M=M, R=R), lambdas)

end

function shakuradelay(lambda::T where T<:Real; M=M, R=R)  # R is accretion rate

    #----------------------------------------------------------------------
    # some constants
    #----------------------------------------------------------------------

    proton_mass = 1.6726231e-24 # in grams
    sun_mass_gm = 1.989e33      # solar mass in grams
    thscs       = 0.665e-24     # cm2 # Thompson scattering cross-section per electron --> (Thompson opacity = thscs / proton_mass = 0.397579 cm2 gram-1)
    g           = 6.67259e-8    # Gravitational constant in cm3 gram-1 s-2
    c           = 2.99792458e10 # speed of light in cm s-1
    year        = 3.15569e7     # seconds in a year
    day         = 86400.0       # seconds in a day
    boltz       = 1.38062e-16   # boltzman constant (erg/K)
    stboltzman  = 5.6696e-5     # Stefan-boltzman constant erg cm-2 s-1 K-4
    H           = 6.6262e-27    # Planck's constant
    angst_cm    = 1e-8          # angstrom in cm
    constant    = (4.0 * pi * g * c * proton_mass)/thscs
    eta         = 0.1           # radiative efficiency or accretion efficiency

    #----------------------------------------------------------------------

    form1 = ((1.0/c) * (2.49 * (boltz/(H*c)) * lambda*angst_cm)^(4.0/3.0))

    L = R * sun_mass_gm / year / 10 * (eta*c^2)


    form2 = (( M * (g*sun_mass_gm/(8.0*pi*stboltzman)) * (10*4*L/(eta*c*c))) ^ (1.0/3.0) ) / day

    tau   = form1 * form2

    return tau

end
