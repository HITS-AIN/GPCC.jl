"""
Uniform prior on delay τ.
It calculates an upper limit on what the delay can be.
This is a theoretical prior based on Eddingtion accretion efficiency.
"""
function uniformpriordelay(; z = z, bhm = bhm, eta = 0.1, edfrac = 10.0)
    
    @printf("Mass is %e\n", bhm)
    @printf("efficiency eta is %f\n", eta)
    @printf("eddington fraction is %f\n", edfrac)

    lum = masslumfunction(;bhm = bhm, edfrac = edfrac, eta = eta)

    f(L,z) = 10.0^(1.559)*(L/10 * 10^(-44))^(0.549) * (1 + z)

    Uniform(0.0, f(lum, z))

end


function masslumfunction(;bhm=1e8,edfrac=10.0,eta=0.1)

    year = 3.15569e7 #seconds in a year
    c    = 2.99792458e10 # speed of light in cm s-1
    g    = 6.67259e-8 #Gravitational constant in cm3 gram-1 s-2
    proton_mass = 1.6726231e-24 # in grams
    thscs       = 0.665e-24 # cm2 # Thompson scattering cross-section per electron --> (Thompson opacity = thscs / proton_mass = 0.397579 cm2 gram-1)
    const_edd   = (4.0 * pi * g * c * proton_mass)/thscs #constant for the eddington luminosity
    sun_mass_gm = 1.989e33 #solar mass in grams
    eddlumin  = const_edd * bhm * sun_mass_gm  #Eddington Luminosity in (erg/s)
    eddrate   = eddlumin/(eta*c^2.0) * year/sun_mass_gm #Eddington mass accretion rate limit
    lum       = edfrac/100.0 * eddlumin

    return lum

end

    
function BLRRatioGivenLum(; z = 0.1, bhm = 1e8, lum = 1e46, eta = 0.1)
        
    #Rin, Rout = BLRRatioGivenLum(z = 0.0258, bhm = 4.6e6, lum = 5.01e43, eta = 0.1)

    proton_mass = 1.6726231e-24 # in grams
    thscs = 0.665e-24 # cm2 # Thompson scattering cross-section per electron --> (Thompson opacity = thscs / proton_mass = 0.397579 cm2 gram-1)
    g = 6.67259e-8 #Gravitational constant in cm3 gram-1 s-2
    c = 2.99792458e10 #speed of light in cm s-1
    year = 3.15569e7 # seconds in a year
    day  = 86400.0 #seconds in a day
    stboltzman = 5.6696e-5 #Stefan-boltzman constant erg cm-2 s-1 K-4
    constan_t = (4.0 * pi * g * c * proton_mass)/thscs
    sun_mass_gm = 1.989e33 #solar mass in grams

    lbol = lum #bolometric luminosity
    mdot = (lbol/(c^2 * eta)) * year/sun_mass_gm #accretion rate in solar masses/year

    Tdust = 1500.0 #Dust temperature in Kelvin.
    bhm = bhm*sun_mass_gm #black hole mass in gram
    mdot = mdot*sun_mass_gm/year  #accretion rate in gram/s

    form1 = 3*g*bhm*mdot
    form2 = 8.0*pi*stboltzman*Tdust^4

    Rin = (form1/form2)^(1.0/3.0) #in cm
    Rin = Rin/c #in seconds
    Rin = Rin/day #in days
    Rin = Rin * (1.0 + z) #delay redshift corrected

    form1 = eta*mdot*c^2
    form2 = 4.0*pi*stboltzman*Tdust^4

    Rout = (form1/form2)^(1.0/2.0) #in cm

    Rout = Rout/c #in seconds
    Rout = Rout/day #in days

    Rout = Rout * (1.0 + z) #delay redshift corrected

    return Rin,Rout
end
    

function BLRRatio(; z = 0.1, bhm = 1e8, edfrac = 10.0, eta = 0.1)
        
    #Rin, Rout = BLRRatio(z = 0.1, bhm = 1e8, edfrac = 10.0, eta = 0.1)

    proton_mass = 1.6726231e-24 # in grams
    thscs = 0.665e-24 # cm2 # Thompson scattering cross-section per electron --> (Thompson opacity = thscs / proton_mass = 0.397579 cm2 gram-1)
    g = 6.67259e-8 #Gravitational constant in cm3 gram-1 s-2
    c = 2.99792458e10 #speed of light in cm s-1
    year = 3.15569e7 # seconds in a year
    day  = 86400.0 #seconds in a day
    stboltzman = 5.6696e-5 #Stefan-boltzman constant erg cm-2 s-1 K-4
    constan_t = (4.0 * pi * g * c * proton_mass)/thscs
    sun_mass_gm = 1.989e33 #solar mass in grams
    eddlumin  = constan_t * bhm * sun_mass_gm
    mdoted = eddlumin/(eta*c^2)
    eddarate = mdoted*year/sun_mass_gm
    mdot = edfrac*eddarate/100.0 #mass accretion rate in solar Masses/year

    Tdust = 1500.0 #Dust temperature in Kelvin.
    bhm = bhm*sun_mass_gm #black hole mass in gram
    mdot = mdot*sun_mass_gm/year  #accretion rate in gram/s

    form1 = 3*g*bhm*mdot
    form2 = 8.0*pi*stboltzman*Tdust^4

    Rin = (form1/form2)^(1.0/3.0) #in cm
    Rin = Rin/c #in seconds
    Rin = Rin/day #in days
    Rin = Rin * (1.0 + z) #delay redshift corrected

    form1 = eta*mdot*c^2
    form2 = 4.0*pi*stboltzman*Tdust^4

    Rout = (form1/form2)^(1.0/2.0) #in cm

    Rout = Rout/c #in seconds
    Rout = Rout/day #in days

    Rout = Rout * (1.0 + z) #delay redshift corrected

    return Rin,Rout
end
    
