"""
This file is to calculate the initial conditions for the jet and the ambient medium that need to be set in PLUTO.

"""
import numpy as np

m_p=1.6726231e-24 #grams
c_light=2.99792458e10 #cm/s
sT=6.65246e-25 #thomson cross section

def gamma(velocity, divided_by_c=False):
    if divided_by_c:
        val=1/np.sqrt(1-velocity**2)
    else:
        val= 1 / np.sqrt(1 - (velocity/c_light) ** 2)
    return val

def calcJetDensPres(lumi_total, gamma_inj, gamma_infinity, opening_angle, inner_jet_radius, unit_dens=1, unit_pres=c_light**2):
    # calculates the density and pressure in real and PLUTO units (same as FLASH units where density
    # is normalized by 1 and pressure is normalized by c**2
    # opening angle should be in degrees
    #for 16TI sim do: calcJetDensPres(5.33e50, 5, 400, 10, 1e9)
    #for 20sp_down ic.calcJetDensPres(2*5.33e50, 5, 400, 10, 1e8) for first period and for FALLBACK_ENVELOPE_20SP_DOWN

    dens=lumi_total/(2*np.pi*c_light**3*inner_jet_radius**2*gamma_inj*gamma_infinity*(1-np.cos(opening_angle*np.pi/180)))

    pres=(0.25*dens*c_light**2)*((gamma_infinity/gamma_inj)-1)

    print('The density is %e in real units and %e in code units (divided by %e g/cm^3)\nThe pressure is %e in real units \
    and %e in code units (divided by %e)\n'%(dens, dens/unit_dens, unit_dens, pres, pres/unit_pres, unit_pres))

    return

def ambientMediumDens(domain_size):
    #domain size must be in cm, prints the min density to have a optical depth of at most 1
    dens=m_p/(sT*domain_size)
    print ('The ambient medium density has to be equal to or smaller than %e to have an optical depth of at \
    most 1\nWith a ratio of 1e-3 between pressure and density times c^2, the pressure needs to be %e'%(dens, dens*c_light**2*1e-3))

    return

def variableLinearEnvelopeJetProfile(tmax=50, jet_on=40,period=10, fraction_constant=0.5,lumi_max=1e50, lumi_init_frac=1, lumi_fraction_decrease=0.05 ):
    #calculates the variabel jet profile where the jet is turned on then off over some period of time. The constant luminosity portion is on for some fraction_constant of the jet period
    #for 20sp_down do ic.variableLinearEnvelopeJetProfile(lumi_max=2*5.33e50, lumi_init_frac=2, period =1, jet_on=20 )

    x=np.linspace(0,tmax, num=10000)
    y=np.zeros_like(x)

    count=0
    period_num=0
    for count in range(x.size):
      #if x[count]> (period_num+1)*period:
        #if (x[count]<jet_on):
        period_num=np.floor(x[count]/period)

        if (x[count]>=period_num*period) & (x[count]<(period_num+fraction_constant)*period) & (x[count]<jet_on):
            factor=lumi_init_frac-period_num*lumi_fraction_decrease
        else:
            if (x[count]<jet_on):
                t=x[count]#-period_num*period
                #print(t, x[count], period_num*period, (fraction_constant)*period)
                factor=(lumi_init_frac-period_num*lumi_fraction_decrease)*(x[count])**(-5/3)/((fraction_constant+period_num)*period)**(-5/3)
            else:
                period_num=np.floor(jet_on-1/period)
                factor = (lumi_init_frac - period_num * lumi_fraction_decrease) * (x[count]) ** (-5 / 3) / (
                          (fraction_constant + period_num) * period) ** (-5 / 3)

        y[count]=lumi_max*factor
      
    return x,y

def variableFallbackEnvelopeJetProfile(tmax=50, jet_on=20, tflat=15, period=1, fraction_constant=0.5, lumi_max=1e50, lumi_init_frac=1):
    #tflat from: https://academic.oup.com/mnras/article/436/2/1867/1155124 and https://iopscience.iop.org/article/10.1086/307790
    #for 20sp_down_global_fallback do: ic.variableFallbackEnvelopeJetProfile(lumi_max=2*5.33e50)

    x=np.linspace(0,tmax, num=1000)
    y=np.zeros_like(x)

    count = 0
    period_num = 0
    for count in range(x.size):
        #if x[count] < jet_on:
        period_num = np.floor(x[count] / period)

        if x[count] <= tflat:
            if (x[count] >= period_num * period) & (x[count] < (period_num + fraction_constant) * period):
                factor = lumi_init_frac
            else:
                factor = (lumi_init_frac)*(x[count]-period_num * period) ** (-5 / 3) / ((fraction_constant) * period) ** (-5 / 3)
        else:
            if (x[count] >= period_num * period) & (x[count] < (period_num + fraction_constant) * period) & (
                    x[count] < jet_on):
                factor = lumi_init_frac*(((period_num + fraction_constant) * period)/tflat)**-(5/3)
            else:
                if (x[count] < jet_on):
                    t = x[count]  # -period_num*period
                    # print(t, x[count], period_num*period, (fraction_constant)*period)
                    factor = (lumi_init_frac*(((period_num + fraction_constant) * period)/tflat)**-(5/3)) * \
                             (x[count]-period_num*period) ** (-5 / 3) / ((fraction_constant) * period) ** (-5 / 3)
                else:
                    cont_value=(lumi_init_frac*(((np.floor(jet_on-1 / period) + fraction_constant) * period)/tflat)**-(5/3)) * \
                             (jet_on-np.floor(jet_on-1 / period)*period) ** (-5 / 3) / ((fraction_constant) * period) ** (-5 / 3) #
                    factor= cont_value*(x[count]/jet_on)**(-5/3)#(lumi_init_frac)*(x[count]/tflat)**-(5/3)*(jet_on/tflat)**-(5/3)*(jet_on-np.floor(jet_on-1 / period)*period) ** (-5 / 3) / ((fraction_constant) * period) ** (-5 / 3)

        y[count] = lumi_max * factor

    return x, y


