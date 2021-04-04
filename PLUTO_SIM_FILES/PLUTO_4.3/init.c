/* ///////////////////////////////////////////////////////////////////// */
/*! 
  \file  
  \brief Relativistic jet propagation.

  This problem sets initial and boundary conditions for a jet
  propagating into a uniform medium with constant density and pressure.
  The computation can be carried using \c CYLINDRICAL or \c SPHERICAL
  coordinates.
  In the first case, the jet enters at the lower z-boundary with speed
  \f$ v_z = \beta\f$ while in the second case, a conical beam with a
  small aperture (\f$\theta = 5^\circ\f$) is injected with the same
  speed from the lower radial boundary.
  At the border of the nozzle, jet values are smoothly joined with
  ambient values using the Profile() function
  (in the current setting the profile is a sharp transition).
  The jet is pressure-matched so that the beam and ambient pressure
  coincide.
 
  The configuration is defined in terms of the following parameters:

  -# <tt>g_inputParam[BETA]</tt>:     the jet velocity;
  -# <tt>g_inputParam[RHO_IN]</tt>:   the jet density;
  -# <tt>g_inputParam[RHO_OUT]</tt>:  the ambient density.
  -# <tt>g_inputParam[PRESS_IN]</tt>: the jet pressure (also equal to ambient
                                     pressure)

  defined in \c pluto.ini.
  The \c TAUB equation of state is used.

  - Configurations #01 and #02 use \c CYLINDRICAL coordinates;
  - Configuration #03 employs \c SPHERICAL coordinates (see snapshot below)

  \image html rhd_jet.03.jpg "Density (log) for configuration #03 at t=200"

  \author A. Mignone (mignone@ph.unito.it)
  \date   Sept 18, 2014
*/
/* ///////////////////////////////////////////////////////////////////// */
#include "pluto.h"


/* ********************************************************************* */
void Init (double *v, double x1, double x2, double x3)
/*
 *  Initalizes the ambient density and pressure to the values given in the pluto.ini file.
 *
 *
 *********************************************************************** */
{
  double scrh;
  #if EOS == IDEAL
   g_gamma = g_inputParam[ADIABATIC_INDEX];
  #endif
  
  EXPAND(v[VX1] = 0.0;, v[VX2] = 0.0;, v[VX3] = 0.0;)//sets the velocities=0 no matter what coordinate system is used
  v[RHO] = g_inputParam[RHO_AMBIENT];
  v[PRS] = g_inputParam[PRESS_AMBIENT];

  //g_smallPressure = v[PRS]/500.0;
}

/* ********************************************************************* */
void InitDomain (Data *d, Grid *grid)
/*! 
 * Assign initial condition by looping over the computational domain.
 * Called after the usual Init() function to assign initial conditions
 * on primitive variables.
 * Value assigned here will overwrite those prescribed during Init().
 * Used to assign the values used for the initial mass distribution
 *
 *********************************************************************** */
{
    int   i, j, k, id;
    double  *x1, *x2, *x3;
    double r;
    
    x1 = grid->xgc[IDIR];  /* -- array pointer to x1 coordinate -- */
    x2 = grid->xgc[JDIR];  /* -- array pointer to x2 coordinate -- */
    x3 = grid->xgc[KDIR];  /* -- array pointer to x3 coordinate -- */

    id=InputDataOpen("STAR_FORMATTED/star_dens_spherical.dbl", "STAR_FORMATTED/star_grid_spherical.out", " ", 0);
    
    TOT_LOOP(k,j,i)
    {
        r=D_EXPAND(x1[i]*x1[i], +x2[j]*x2[j], +x3[k]*x3[k]);
        r=sqrt(r);
        
        //if r is less than the radius of the star do the interpolation
        if (r<=g_inputParam[R_STAR])
        {
            d->Vc[RHO][k][j][i]=InputDataInterpolate(id, x1[i], x2[j], x3[k]);
        }
    }
    InputDataClose(id);
    
    id=InputDataOpen("STAR_FORMATTED/star_pres_spherical.dbl", "STAR_FORMATTED/star_grid_spherical.out", " ", 0);
    
    TOT_LOOP(k,j,i)
    {
        r=D_EXPAND(x1[i]*x1[i], +x2[j]*x2[j], +x3[k]*x3[k]);
        r=sqrt(r);
        
        //if r is less than the radius of the star do the interpolation
        if (r<=g_inputParam[R_STAR])
        {
            d->Vc[PRS][k][j][i]=InputDataInterpolate(id, x1[i], x2[j], x3[k]);
        }
    }
    InputDataClose(id);
    
    //assign velocities to be 0 everywhere
    TOT_LOOP(k,j,i)
    {
        EXPAND(d->Vc[VX1][k][j][i] = 0.0;, d->Vc[VX2][k][j][i] = 0.0;, d->Vc[VX3][k][j][i] = 0.0;)
    }
    
}

/* ********************************************************************* */
void Analysis (const Data *d, Grid *grid)
/* 
 *
 *
 *********************************************************************** */
{

}
/* ********************************************************************* */
void UserDefBoundary (const Data *d, RBox *box, int side, Grid *grid) 
/*
 *
 *
 *********************************************************************** */
{
    
  int   i, j, k, nv;
  double  *x1, *x2, *x3;
  double r, rmin, vj[NVAR];
    double angle, factor;
    double period_num, time;
    double cont_value;
  //double cs, Tj, vj[NVAR];
  //double  r, vjet[NVAR], vout[NVAR];

  x1 = grid->xgc[IDIR];  /* -- array pointer to x1 coordinate -- */
  x2 = grid->xgc[JDIR];  /* -- array pointer to x2 coordinate -- */
  x3 = grid->xgc[KDIR];  /* -- array pointer to x3 coordinate -- */
    
// this is a nother way from another example file
    //if the time is less than some value the jet is on
    /*
    if (g_time<g_inputParam[JET_ON]+1)
    {
        if (side == X1_BEG)
        {     // -- select the boundary side --
            BOX_LOOP(box,k,j,i)
            {   // -- Loop over boundary zones --
                //angle=atan(x1[i]/(x2[j]));
                if ( (x2[j]<=g_inputParam[THETA_JET]*CONST_PI/180))
                {   // -- set jet values for x within opening angle of jet
                    d->Vc[RHO][k][j][i] = g_inputParam[RHO_JET];
                    d->Vc[VX1][k][j][i] = g_inputParam[VEL_JET];
                    d->Vc[VX2][k][j][i] = 0.0;
                    //#if GEOMETRY == SPHERICAL
                    //    d->Vc[VX3][k][j][i] = vx3;
                    //#endif
                    d->Vc[PRS][k][j][i] = g_inputParam[PRESS_JET];
                }
                else
                {                // -- reflective boundary for r > 1 --
                    VAR_LOOP(nv) d->Vc[nv][k][j][i] = d->Vc[nv][k][j][2*IBEG-i-1];
                    //d->Vc[VX1][k][j][i] *= -1.0;
                    d->Vc[VX1][k][j][i] *= -1.0;
                    //#if GEOMETRY == SPHERICAL
                    //    d->Vc[VX3][k][j][i] *= -1.0;
                    //#endif
                }
            }
        }
    }
    else
    {
        //if the time is past the jet "on" time, the jet is off
        if (side == X1_BEG)
        {     // -- select the boundary side --
            BOX_LOOP(box,k,j,i)
            {   // -- Loop over boundary zones --
                if ( (x2[j]<=g_inputParam[THETA_JET]*CONST_PI/180))
                {   // -- set jet values for x within opening angle of jet
                    //for the jet being shut off
                    d->Vc[RHO][k][j][i] = g_inputParam[RHO_JET]*pow(g_time-g_inputParam[JET_ON], -5.0/3.0);
                    d->Vc[VX1][k][j][i] = g_inputParam[VEL_JET]; //*pow(g_time-100, -5.0/3.0);
                    d->Vc[VX2][k][j][i] = 0.0;
                    //#if GEOMETRY == SPHERICAL
                    //    d->Vc[VX3][k][j][i] = vx3;
                    //#endif
                    d->Vc[PRS][k][j][i] = g_inputParam[PRESS_JET]*pow(g_time-g_inputParam[JET_ON], -5.0/3.0);
                    
                    printf("Jet BCs at end: %e, %e, %e, %e, %e\n", g_time, pow(g_time-g_inputParam[JET_ON], -5.0/3.0),d->Vc[RHO][k][j][i] ,d->Vc[VX1][k][j][i], d->Vc[PRS][k][j][i])
                }
                else
                {                // -- reflective boundary everywhere --
                    VAR_LOOP(nv) d->Vc[nv][k][j][i] = d->Vc[nv][k][j][2*IBEG-i-1];
                    //d->Vc[VX1][k][j][i] *= -1.0;
                    d->Vc[VX1][k][j][i] *= -1.0;
                    //#if GEOMETRY == SPHERICAL
                    //    d->Vc[VX3][k][j][i] *= -1.0;
                    //#endif
                }
            }
        }
    }
    */
    
    //to define the min density/pressure
    if (side == 0)
    {
        TOT_LOOP(k,j,i)
        {
           if ((d->Vc[RHO][k][j][i] < g_inputParam[RHO_AMBIENT]) || (d->Vc[PRS][k][j][i] < g_inputParam[PRESS_AMBIENT]))
           {
              d->Vc[RHO][k][j][i] = g_inputParam[RHO_AMBIENT];
              d->Vc[PRS][k][j][i] < g_inputParam[PRESS_AMBIENT];
           }
        }
    }
    
    //for  constant injected jet
    if (g_time<g_inputParam[JET_ON]+1)
    {
        factor=1;
    }
    else
    {
        factor=pow(g_time-g_inputParam[JET_ON], -5.0/3.0);//not quite t^-5/3, but gets there asymptotically
    }
    //for  constant injected jet

    
    //for a variable injected jet with linealy decreasing envelope, can also be for constant jet if JET_FRAC_ON is set to 1
    period_num=floor(g_time/g_inputParam[JET_PERIOD]);

    if ((g_time>=period_num*g_inputParam[JET_PERIOD]) && (g_time<(period_num+g_inputParam[JET_FRAC_ON])*g_inputParam[JET_PERIOD] ) && (g_time<g_inputParam[JET_ON]))
    {
        factor=g_inputParam[JET_LUMI_INIT_FRAC]-period_num*g_inputParam[JET_LUMI_FRAC_DEC];
    }
    else
    {
        if (g_time<g_inputParam[JET_ON])
        {
                time=g_time; //-period_num*g_inputParam[JET_PERIOD];
                factor=(g_inputParam[JET_LUMI_INIT_FRAC]-period_num*g_inputParam[JET_LUMI_FRAC_DEC])*pow(time, -5.0/3.0)/pow((g_inputParam[JET_FRAC_ON]+period_num)*g_inputParam[JET_PERIOD], -5.0/3.0);
        }
        else
        {
            period_num=floor((g_inputParam[JET_ON]-1)/g_inputParam[JET_PERIOD]);
            factor=(g_inputParam[JET_LUMI_INIT_FRAC]-period_num*g_inputParam[JET_LUMI_FRAC_DEC])*pow(g_time, -5.0/3.0)/pow((g_inputParam[JET_FRAC_ON]+period_num)*g_inputParam[JET_PERIOD], -5.0/3.0);

        }
    }
    //for a variable injected jet with linealy decreasing envelope, can also be for constant jet if JET_FRAC_ON is set to 1


    //for a variable injected jet with a flat global envelope followed by a t^-5/3 envelope
    period_num=floor(g_time/g_inputParam[JET_PERIOD]);
    
    if (g_time<=g_inputParam[JET_TFLAT])
    {
        if ((g_time>=period_num*g_inputParam[JET_PERIOD]) && (g_time<(period_num+g_inputParam[JET_FRAC_ON])*g_inputParam[JET_PERIOD] ))
        {
            factor=g_inputParam[JET_LUMI_INIT_FRAC];
        }
        else
        {
            factor=g_inputParam[JET_LUMI_INIT_FRAC]*pow(g_time-period_num*g_inputParam[JET_PERIOD], -5.0/3.0)/pow(g_inputParam[JET_FRAC_ON]*g_inputParam[JET_PERIOD], -5.0/3.0);
        }
    }
    else
    {
        if ((g_time>=period_num*g_inputParam[JET_PERIOD]) && (g_time<(period_num+g_inputParam[JET_FRAC_ON])*g_inputParam[JET_PERIOD] ) && (g_time<g_inputParam[JET_ON]))
        {
            factor=g_inputParam[JET_LUMI_INIT_FRAC]*pow((period_num+g_inputParam[JET_FRAC_ON])*g_inputParam[JET_PERIOD]/g_inputParam[JET_TFLAT], -5.0/3.0);
        }
        else
        {
            if (g_time<g_inputParam[JET_ON])
            {
                factor = g_inputParam[JET_LUMI_INIT_FRAC]*pow((period_num+g_inputParam[JET_FRAC_ON])*g_inputParam[JET_PERIOD]/g_inputParam[JET_TFLAT], -5.0/3.0)*pow(g_time-period_num*g_inputParam[JET_PERIOD], -5.0/3.0)/pow(g_inputParam[JET_FRAC_ON]*g_inputParam[JET_PERIOD], -5.0/3.0);
                
            }
            else
            {
                period_num=floor((g_inputParam[JET_ON]-1)/g_inputParam[JET_PERIOD]);
                cont_value=g_inputParam[JET_LUMI_INIT_FRAC]*pow((period_num+g_inputParam[JET_FRAC_ON])*g_inputParam[JET_PERIOD]/g_inputParam[JET_TFLAT], -5.0/3.0)*pow(g_inputParam[JET_ON]-period_num*g_inputParam[JET_PERIOD], -5.0/3.0)/pow(g_inputParam[JET_FRAC_ON]*g_inputParam[JET_PERIOD], -5.0/3.0);
                
                factor= cont_value*pow(g_time/g_inputParam[JET_ON], -5.0/3.0);

            }
        }
    }
    //for a variable injected jet with a flat global envelope followed by a t^-5/3 envelope

    
    
    if (side == X1_BEG)
    {     // -- select the boundary side --
        BOX_LOOP(box,k,j,i)
        {   // -- Loop over boundary zones --
            //angle=atan(x1[i]/(x2[j]));
            if ( (x2[j]<=g_inputParam[THETA_JET]*CONST_PI/180))
            {   // -- set jet values for x within opening angle of jet
                d->Vc[RHO][k][j][i] = g_inputParam[RHO_JET]*factor;
                d->Vc[VX1][k][j][i] = g_inputParam[VEL_JET]; //shoild this be *factor? NO it shouldnt
                d->Vc[VX2][k][j][i] = 0.0;
                //#if GEOMETRY == SPHERICAL
                //    d->Vc[VX3][k][j][i] = vx3;
                //#endif
                d->Vc[PRS][k][j][i] = g_inputParam[PRESS_JET]*factor;
            }
            else
            {                // -- reflective boundary for r > 1 --
                VAR_LOOP(nv) d->Vc[nv][k][j][i] = d->Vc[nv][k][j][2*IBEG-i-1];
                //d->Vc[VX1][k][j][i] *= -1.0;
                d->Vc[VX1][k][j][i] *= -1.0;
                //#if GEOMETRY == SPHERICAL
                //    d->Vc[VX3][k][j][i] *= -1.0;
                //#endif
            }
        }
    }
    
    
    
    
}


