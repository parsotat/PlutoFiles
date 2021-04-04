#include <cstdio>
#include <string>
using std::string;

#include "PatchPluto.H"
#include "LoHiSide.H"

extern "C" {
#include "pluto.h"
}

static void computeRefVar(double ***UU[], double ***q, double, RBox *Ubox);

#if (EOS != ISOTHERMAL) && (ENTROPY_SWITCH == NO)
 #ifndef CHOMBO_REF_VAR  
  #define CHOMBO_REF_VAR ENG 
 #endif
#else
 #ifndef CHOMBO_REF_VAR  
  #define CHOMBO_REF_VAR RHO
 #endif
#endif

#define REF_CRIT 2   /* 1 == first derivative, 2 == second derivative */

/* ************************************************************************* */
void PatchPluto::computeRefGradient(FArrayBox& gFab, FArrayBox& UFab, 
                                    const FArrayBox& a_dV, const Box& b)
/*!
 * Tag zones for refinement using gradient of the conservative 
 * variables.
 * The gradient is computed by standard finite differences using
 *
 * - REF_CRIT equal to 1 --> compute (normalized) gradient using 1st 
 *                           derivative of the solution;
 * - REF_CRIT equal to 2 --> compute (normalized) gradient using 2nd 
 *                           derivative of the solution (default);
 *                           This approach is based on Lohner (1987).
 *
 * Zones will be flagged for refinement whenever grad[k][j][i] exceeds 
 * the threshold value specified by the 'Refine_thresh' parameter read in
 * pluto.ini.
 *
 * Derivatives are computed using the conserved variable
 * U[CHOMBO_REF_VAR] 
 * where CHOMBO_REF_VAR is taken to be energy density (default).
 * However, by setting CHOMBO_REF_VAR = -1, you can provide your own 
 * physical variable through the function computeRefVar().
 * 
 * \authors C. Zanni   (zanni@oato.inaf.it)\n
 *          A. Mignone (mignone@ph.unito.it)
 * \date    Oct 11, 2012
 *************************************************************************** */
{
  CH_assert(m_isDefined);

  int nv, i, j, k;
  double x1, dqx_p, dqx_m, dqx, d2qx, den_x;
  double x2, dqy_p, dqy_m, dqy, d2qy, den_y;
  double x3, dqz_p, dqz_m, dqz, d2qz, den_z;
  double gr1, gr2, eps = 0.01;
#if CHOMBO_REF_VAR == -1
  double ***UU[NVAR];
#endif
  double ***q, ***grad;
  RBox  Ubox, Gbox;
    
    double argl, argr, x_l, x_r, x, dx, xgc, r, r_jet_back, r_jet_front, r_jet_avg, delta_theta_ref;
    int max_level, count_level, count_level_angle;
/* -- check ref criterion -- */

#if REF_CRIT != 1 && REF_CRIT != 2
  print ("! TagCells.cpp: Refinement criterion not valid\n");
  QUIT_PLUTO(1);
#endif

/* -----------------------------------------------------
   1. The solution array U is defined on the box 
      [Uib, Uie] x [Ujb, Uje] x [Ukb, Uke], which 
      differs from that of gFab ([Gib,...Gke]), 
      typically one point larger in each direction. 
   ----------------------------------------------------- */
    
  Ubox.jbeg = Ubox.jend = Ubox.kbeg = Ubox.kend = 0;
  Gbox.jbeg = Gbox.jend = Gbox.kbeg = Gbox.kend = 0;

  D_EXPAND(Ubox.ibeg = UFab.loVect()[IDIR]; Ubox.iend = UFab.hiVect()[IDIR]; ,
           Ubox.jbeg = UFab.loVect()[JDIR]; Ubox.jend = UFab.hiVect()[JDIR]; ,
           Ubox.kbeg = UFab.loVect()[KDIR]; Ubox.kend = UFab.hiVect()[KDIR]; );

  D_EXPAND(Gbox.ibeg = gFab.loVect()[IDIR]; Gbox.iend = gFab.hiVect()[IDIR]; ,
           Gbox.jbeg = gFab.loVect()[JDIR]; Gbox.jend = gFab.hiVect()[JDIR]; ,
           Gbox.kbeg = gFab.loVect()[KDIR]; Gbox.kend = gFab.hiVect()[KDIR]; );

/* --------------------------------------------------------
   2. Input solution array (UFab.dataPtr(nv)) is defined 
      as dV*U/dx^3, where U is an array of conservative 
      variables and dV is the zone volume. 
      To obtain U we must divide by volume.
   -------------------------------------------------------- */

#if CHOMBO_REF_VAR == -1
  FArrayBox tmpU(UFab.box(),NVAR);
  tmpU.copy(UFab);
#else
  FArrayBox tmpU(UFab.box(),1);
  tmpU.copy(UFab,CHOMBO_REF_VAR,0);
#endif 

#if GEOMETRY != CARTESIAN

  #if CHOMBO_REF_VAR == -1

    for (nv = 0; nv < NVAR; nv++) tmpU.divide(a_dV,0,nv);
    #if CHOMBO_CONS_AM == YES
      #if ROTATING_FRAME == YES
        Box curBox = UFab.box();
        for(BoxIterator bit(curBox); bit.ok(); ++bit) {
          const IntVect& iv = bit();
          tmpU(iv,iMPHI) /= a_dV(iv,1);
          tmpU(iv,iMPHI) -= tmpU(iv,RHO)*a_dV(iv,1)*g_OmegaZ;
        }
      #else
        tmpU.divide(a_dV,1,iMPHI);
      #endif
    #endif

  #else

    tmpU.divide(a_dV,0,0);

  #endif

#else

   if (g_stretch_fact != 1.) tmpU /= g_stretch_fact;

#endif // GEOMETRY == CARTESIAN

/* ---------------------------------------------
   3. Set refinement variable
   --------------------------------------------- */

#if CHOMBO_REF_VAR == -1
  for (nv = 0; nv < NVAR; nv++)
    UU[nv] = ArrayBoxMap(Ubox.kbeg, Ubox.kend,
                         Ubox.jbeg, Ubox.jend,
                         Ubox.ibeg, Ubox.iend, tmpU.dataPtr(nv));
  q = ArrayBox(Ubox.kbeg, Ubox.kend, Ubox.jbeg, Ubox.jend, Ubox.ibeg, Ubox.iend);
  computeRefVar(UU, q, m_dx, &Ubox);
#else
  q = ArrayBoxMap(Ubox.kbeg, Ubox.kend,
                  Ubox.jbeg, Ubox.jend,
                  Ubox.ibeg, Ubox.iend, tmpU.dataPtr(0));
#endif
  
  grad = ArrayBoxMap(Gbox.kbeg, Gbox.kend, 
                     Gbox.jbeg, Gbox.jend, 
                     Gbox.ibeg, Gbox.iend, gFab.dataPtr(0));

/* ----------------------------------------------------------------
   4. Main spatial loop for zone tagging based on 1st 
     (REF_CRIT = 1) or 2nd (REF_CRIT = 2) derivative error norm. 
   ---------------------------------------------------------------- */

  BOX_LOOP(&Gbox, k, j, i){
    x3 = (k + 0.5)*m_dx*g_x3stretch + g_domBeg[KDIR];
    x2 = (j + 0.5)*m_dx*g_x2stretch + g_domBeg[JDIR];
#if CHOMBO_LOGR == NO
    x1 = (i + 0.5)*m_dx          + g_domBeg[IDIR];
#else
    double xl = g_domBeg[IDIR] + i*m_dx;
    double xr = xl + m_dx;
    x1 = g_domBeg[IDIR]*0.5*(exp(xr)+exp(xl));
#endif
      
    

    D_EXPAND(dqx_p =    q[k][j][i+1] - q[k][j][i];
             dqx_m = - (q[k][j][i-1] - q[k][j][i]);  ,
             dqy_p =    q[k][j+1][i] - q[k][j][i];
             dqy_m = - (q[k][j-1][i] - q[k][j][i]);  ,
             dqz_p =    q[k+1][j][i] - q[k][j][i];
             dqz_m = - (q[k-1][j][i] - q[k][j][i]);)

  /* --------------------------------------------------------------
      Physical boundary values are not up to date and should be 
      excluded from gradient computation. 
      In this case, left and right derivatives are set equal to 
      each other. This will not trigger refinement in the leftmost 
      and rightmost internal zones (using 2nd derivative) but we 
      really don't care since buffer size will do the job.
     -------------------------------------------------------------- */
      
    D_EXPAND(if (i == 0) dqx_m = dqx_p;  ,
             if (j == 0) dqy_m = dqy_p;  ,
             if (k == 0) dqz_m = dqz_p;)

    D_EXPAND(if (i == m_domain.size(IDIR)-1) dqx_p = dqx_m;  ,
             if (j == m_domain.size(JDIR)-1) dqy_p = dqy_m;  ,
             if (k == m_domain.size(KDIR)-1) dqz_p = dqz_m;)

  /* -----------------------------------------------
         Compute gradient using 1st derivative 
      ---------------------------------------------- */

#if REF_CRIT == 1
    D_EXPAND(dqx = dqx_p + dqx_m;  ,
             dqy = dqy_p + dqy_m;  ,
             dqz = dqz_p + dqz_m;)

    D_EXPAND(den_x = fabs(q[k][j][i+1]) + fabs(q[k][j][i-1]);  ,
             den_y = fabs(q[k][j+1][i]) + fabs(q[k][j-1][i]);  ,
             den_z = fabs(q[k+1][j][i]) + fabs(q[k-1][j][i]);)

    gr1  = D_EXPAND(dqx*dqx, + dqy*dqy, + dqz*dqz);
    gr1 /= D_EXPAND(den_x*den_x, + den_y*den_y, + den_z*den_z);

    grad[k][j][i] = sqrt(gr1);
#endif

  /* -----------------------------------------------
         Compute gradient using 2nd derivative 
      ---------------------------------------------- */

#if REF_CRIT == 2
    D_EXPAND(d2qx = dqx_p - dqx_m;  ,
             d2qy = dqy_p - dqy_m;  ,
             d2qz = dqz_p - dqz_m;)

    D_EXPAND(
      den_x = 2.0*fabs(q[k][j][i]) + fabs(q[k][j][i+1]) + fabs(q[k][j][i-1]);
      den_x = fabs(dqx_p) + fabs(dqx_m) + eps*den_x;    ,

      den_y = 2.0*fabs(q[k][j][i]) + fabs(q[k][j+1][i]) + fabs(q[k][j-1][i]);
      den_y = fabs(dqy_p) + fabs(dqy_m) + eps*den_y;    ,

      den_z = 2.0*fabs(q[k][j][i]) + fabs(q[k+1][j][i]) + fabs(q[k-1][j][i]);
      den_z = fabs(dqz_p) + fabs(dqz_m) + eps*den_z;
    )

    gr2  = D_EXPAND(d2qx*d2qx,   + d2qy*d2qy,   + d2qz*d2qz);
    gr2 /= D_EXPAND(den_x*den_x, + den_y*den_y, + den_z*den_z);

    grad[k][j][i] = sqrt(gr2);
#endif
    
      //try to do time dependent refinement, based on where jet should be located
      //calculate radius for logrithmic radial grid (from pyPLUTO), x2 (theta) is already calculated above
      r=g_domBeg[IDIR]*0.5*(exp(i*m_dx)+exp((i+1)*m_dx )  );


      //set CHOMBO radial refinement, can get tricky when moving from one evel of refinement to the next, so they may need to be overlapped once the jet gets to the region o space that you want to be a higher refinement
    /*
    if (r<5 && r>0 &&  m_level>1)
    {
    grad[k][j][i] = 0.0000;
    }
    else if (r<50 && r>=5 && m_level>2)
    {
    grad[k][j][i] = 0.0000;
    }
    else if (r<500 && r>=50 && m_level>2)
    {
    grad[k][j][i] = 0.0000;
    }
    else if (r<5000 && r>=500 && m_level>3)
    {
    grad[k][j][i] = 0.0000;
    }
    else if (r<550000 && r>=5000 && m_level>4)
    {
    grad[k][j][i] = 0.0000;
    }


    if (x2>20*3.1415/180  && m_level>2) //set maximum CHOMBO refinement in theta
    {
     grad[k][j][i] = 0.0000; //if outside of the 20 degree cone, dont refine past a level of 3
    }
    */
      
      //this is a better way of adaptively refining cells following the jet propagation!!!!!!!!!!!!!!!!!
      //front and back of jet at
      r_jet_front=CONST_c*g_time*UNIT_LENGTH/UNIT_VELOCITY/UNIT_LENGTH;
      r_jet_back=r_jet_front-1.5*CONST_c*g_inputParam[JET_ON]*UNIT_LENGTH/UNIT_VELOCITY/UNIT_LENGTH;
      
      if (r_jet_front<=g_inputParam[R_STAR]) //if the jet hasnt broken out of the star yet, make sure that the star is fully resolved
      {
          r_jet_front=g_inputParam[R_STAR];
          r_jet_back=0;
      }
      
      r_jet_avg=0.5*(r_jet_front+r_jet_back);
      
      //printf("front %e, back %e, avg %e\n", r_jet_front, r_jet_back, r_jet_avg);
      
      if (r_jet_avg<80)
      {
          max_level=3;
      }
      else if (r_jet_avg<500 && r_jet_avg>=80)
      {
          max_level=4;
      }
      else if (r_jet_avg<5000 && r_jet_avg>=500)
      {
          max_level=5;
      }
      else if (r_jet_avg<50000 && r_jet_avg>=5000)
      {
          max_level=6;
      }
      
      count_level=max_level;
      while (!(r<(1+0.05*max_level-0.05*count_level)*r_jet_front && r>(1-0.05*max_level+0.05*count_level)*r_jet_back))
      {
          count_level--; //identify the refinement level that the block should have based on the limits for each region of refinement
      }
            
      count_level_angle=max_level;
      delta_theta_ref=10*3.14/180;
      while (!(x2>(max_level-count_level_angle)*delta_theta_ref && x2<(1+max_level-count_level_angle)*delta_theta_ref  ))
      {
          count_level_angle--;
      }
      
      if (r_jet_back<=g_inputParam[R_STAR] && r<=g_inputParam[R_STAR])
      {
          count_level_angle=3;//if the jet is still in the star make sure that the whole star is resolved to the max with no change in refinement in theta
      }
      
      if (g_time==0)
      {
          count_level=3;//do this as the default when PLUTO starts up again, otherwise get weird behavior
          count_level_angle=3;
      }
      
      
      if (m_level>count_level-1 || m_level>count_level_angle-1)
      {
          grad[k][j][i] = 0.0000;
      }

  }

/* --------------------------------------------------------------
   6. Free array
   -------------------------------------------------------------- */
   
  FreeArrayBoxMap(grad, Gbox.kbeg, Gbox.kend,
                        Gbox.jbeg, Gbox.jend,
                        Gbox.ibeg, Gbox.iend);

#if CHOMBO_REF_VAR == -1
  for (nv = 0; nv < NVAR; nv++){
    FreeArrayBoxMap(UU[nv], Ubox.kbeg, Ubox.kend,
                            Ubox.jbeg, Ubox.jend,
                            Ubox.ibeg, Ubox.iend);
  }
  FreeArrayBox(q, Ubox.kbeg, Ubox.jbeg, Ubox.ibeg);
#else
  FreeArrayBoxMap(q, Ubox.kbeg, Ubox.kend,
                     Ubox.jbeg, Ubox.jend,
                     Ubox.ibeg, Ubox.iend);
#endif

}

/* ********************************************************************* */
void computeRefVar(double ***UU[], double ***q, double dx, RBox *Ubox)
/*!
 * Compute a user-defined array q(U) function of the conserved
 * variables.
 *
 *
 *********************************************************************** */
{
  int nv, i, j, k;
  double us[NVAR], vs[NVAR];

  BOX_LOOP(Ubox, k, j, i) {
    VAR_LOOP(nv) us[nv] = UU[nv][k][j][i]; 
    q[k][j][i] = us[RHO];
  }

}
#undef REF_CRIT
