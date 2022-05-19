/*==========================================================
 * F_comp_2d.c - grid to particle interpolation
 *
 * Computes deformation gradient F 
 *
 * The calling syntax is:
 *
 *		[Fup1,Fup2,Fup3,Fup4] = F_comp_2d(Ox,Oy, hx,hy, nx,ny, Xp,Yp, dt, vx,vy, F1,F2,F3,F4)
 *
 * This is a MEX-file for MATLAB.
 *
 *========================================================*/

#include "mex.h"
#include <stdio.h>
#include <math.h>

// Weight Function
double N_weight(double x) 
{
    x = fabs(x);     // *** abs function is for ints only (in C) ***
    double N = 0.0;
    if (x<0.5)
        N += 0.75 - x*x;
    else if (x<1.5)
        N += 0.5*(1.5-x)*(1.5-x);

    return N;
}

double Nx_weight(double x) 
{
   double Nx = 0.0;
   if (fabs(x)<0.5) // *** abs function is for ints only (in C) ***
      Nx += -2*x;
   else if (x<1.5 && x>=0.5)
      Nx += -1.5 +x;
   else if (x>-1.5 && x<=-0.5)
      Nx += 1.5 +x;   
   return Nx;
}

// Deformation gradient computation routine
void F_comp(const double Ox, const double Oy,                       // scalars
            const double hx, const double hy, 
            size_t nx, size_t ny, size_t NP,
            double *Xp, double *Yp,                                 // input
            double dt, double *vx, double *vy,
            double *F1, double *F2, double *F3, double *F4,  
            double *Fup1, double *Fup2, double *Fup3, double *Fup4  // output       
            )
{
    // Constants
    const int d = 2;
   
    // Temp vars
    double fp_curr1, fp_curr2, fp_curr3, fp_curr4; 
    double xp_curr, yp_curr;
    double xg_curr, yg_curr, dx, dy;
    int I, J;
    int l, k, gidx;
    double Nx, Ny;
   
    // EVERYTHING CONDENSED IN THIS LOOP OVER PARTICLES
    for (int pidx=0; pidx<NP; pidx++) {
        fp_curr1 = 0.0; 
        fp_curr2 = 0.0;
        fp_curr3 = 0.0; 
        fp_curr4 = 0.0;
        xp_curr = Xp[pidx];
        yp_curr = Yp[pidx];

        // Index of reference grid node for particle p 
        I = floor( (xp_curr-Ox) / hx );  
        J = floor( (yp_curr-Oy) / hy );  
      
        // 1st loop over stencil nodes
        for (int j=0; j<49; j++) {
            k = floor(j/7);      // k: 0 to 6
            l = j - (7*k);       // l: 0 to 6
            k = I + k - 3;       // k: Ip + (-3 to 3)
            l = J + l - 3;       // l: Jp + (-3 to 3)
         
            // Check if stencil indices are valid 
            if (k>=0 && k<nx+1 && l>=0 && l<ny+1) {
                // Compute gidx of (k,l) 
                gidx = (k)*(ny+1) + l; 
                xg_curr = Ox + hx*k;
                yg_curr = Oy + hy*l;
            
                // Compute Interpolation Weight
                dx = (xp_curr - xg_curr)/hx; 
                dy = (yp_curr - yg_curr)/hy; 
            
            // PIC transfer
            Nx = Nx_weight(dx) * N_weight(dy) / hx;
            Ny = N_weight(dx) * Nx_weight(dy) / hy;
            fp_curr1 += vx[gidx]*Nx;
            fp_curr2 += vy[gidx]*Nx;
            fp_curr3 += vx[gidx]*Ny;
            fp_curr4 += vy[gidx]*Ny;
            }
        }
        // Ft = Id + dt*Ft
        fp_curr1 = dt*fp_curr1 + 1;       
        fp_curr2 = dt*fp_curr2;
        fp_curr3 = dt*fp_curr3;
        fp_curr4 = dt*fp_curr4 + 1;
      
        // Fup = Ft * Fkp
        Fup1[pidx] = fp_curr1*F1[pidx] + fp_curr3*F2[pidx];
        Fup2[pidx] = fp_curr2*F1[pidx] + fp_curr4*F2[pidx];
        Fup3[pidx] = fp_curr1*F3[pidx] + fp_curr3*F4[pidx];
        Fup4[pidx] = fp_curr2*F3[pidx] + fp_curr4*F4[pidx];
    }
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    // Declare variables
    double *Xp, *Yp;                    // Particle positions (NP x 1)
    double dt;                          // time step
    double *vx, *vy;                    // grid vel v=(u+w)   (NG x 1)
    double *F1, *F2, *F3, *F4;          // def. grad. Fpk     (NP x 1)
    double *Fup1, *Fup2, *Fup3, *Fup4;  // F update           (NP x 1)
    //const double Ox; const double Oy;   // "origin" of Eul. grid
    //const double hx; const double hy;   // grid spacing
    size_t nx; size_t ny;               // Eul. grid partitions
    size_t NP;                          // # particles

    /* check for proper number of arguments */
    if (nrhs!=15) {
        mexErrMsgIdAndTxt("F_comp_2d:nrhs","15 inputs required.");
    }
    if (nlhs!=4) {
        mexErrMsgIdAndTxt("F_comp_2d:nlhs","Four outputs required.");
    }

    /* make sure the 6th input argument is type double */
    if (!mxIsDouble(prhs[6]) ||  mxIsComplex(prhs[6])) {
        mexErrMsgIdAndTxt("F_comp_2d:notDouble","Xp must be type double.");
    }

    /* check that number of cols in input arguments is 1 */
    if (mxGetN(prhs[6])!=1 || mxGetN(prhs[7])!=1 ) {
        mexErrMsgIdAndTxt("F_comp_2d:notRowVector","Xp, Yp must be a column vectors.");
    }
   
    /* check that number of rows in input arguments is 1 */
    if (mxGetM(prhs[9])!=1 || mxGetM(prhs[10])!=1 || 
        mxGetM(prhs[11])!=1 || mxGetM(prhs[12])!=1 || 
        mxGetM(prhs[13])!=1 || mxGetM(prhs[14])!=1 ) {
        mexErrMsgIdAndTxt("F_comp_2d:notRowVector","vx, vy, F1, F2, F3, F4 must be a row vectors.");
    }
   
    /* make sure the following input arguments are scalar */
    if (mxGetNumberOfElements(prhs[0])!=1 || mxGetNumberOfElements(prhs[1])!=1 ||
        mxGetNumberOfElements(prhs[2])!=1 || mxGetNumberOfElements(prhs[3])!=1 ||
        mxGetNumberOfElements(prhs[4])!=1 || mxGetNumberOfElements(prhs[5])!=1 ||
        mxGetNumberOfElements(prhs[8])!=1) {
        mexErrMsgIdAndTxt("F_comp_2d:notScalar","Inputs Ox, Oy, hx, hy, nx, ny, dt must be scalars.");
    }
   
    /* get the values of the scalar inputs  */
    const double Ox = mxGetScalar(prhs[0]);
    const double Oy = mxGetScalar(prhs[1]);
    const double hx = mxGetScalar(prhs[2]);
    const double hy = mxGetScalar(prhs[3]);
    nx = mxGetScalar(prhs[4]);
    ny = mxGetScalar(prhs[5]);
    dt = mxGetScalar(prhs[8]);
   
    /* create a pointer to the real data in Xp, Yp, vx, vy, F1, F2, F3, F4 */
    #if MX_HAS_INTERLEAVED_COMPLEX  // get a C error if we don't do this if/else statement
    Xp = mxGetDoubles(prhs[6]);
    Yp = mxGetDoubles(prhs[7]);
    vx = mxGetDoubles(prhs[9]);
    vy = mxGetDoubles(prhs[10]);
    F1 = mxGetDoubles(prhs[11]);
    F2 = mxGetDoubles(prhs[12]);
    F3 = mxGetDoubles(prhs[13]);
    F4 = mxGetDoubles(prhs[14]);
    #else
    Xp = mxGetPr(prhs[6]);          // matlab recommends not to use mxGetPr
    Yp = mxGetPr(prhs[7]);
    vx = mxGetPr(prhs[9]);
    vy = mxGetPr(prhs[10]);
    F1 = mxGetPr(prhs[11]);
    F2 = mxGetPr(prhs[12]);
    F3 = mxGetPr(prhs[13]);
    F4 = mxGetPr(prhs[14]);
    #endif
   
    /* get dimensions of the input matrix */
    NP = mxGetM(prhs[6]);

    /* create the output matrices */
    plhs[0] = mxCreateDoubleMatrix(1,(mwSize)NP,mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1,(mwSize)NP,mxREAL);
    plhs[2] = mxCreateDoubleMatrix(1,(mwSize)NP,mxREAL);
    plhs[3] = mxCreateDoubleMatrix(1,(mwSize)NP,mxREAL);

    /* get a pointer to the real data in the output matrix */
    #if MX_HAS_INTERLEAVED_COMPLEX
    Fup1 = mxGetDoubles(plhs[0]);
    Fup2 = mxGetDoubles(plhs[1]);
    Fup3 = mxGetDoubles(plhs[2]);
    Fup4 = mxGetDoubles(plhs[3]);
    #else
    Fup1 = mxGetPr(plhs[0]); 
    Fup2 = mxGetPr(plhs[1]);
    Fup3 = mxGetPr(plhs[2]); 
    Fup4 = mxGetPr(plhs[3]);
    #endif
   
    /* call the computational routine */
    F_comp(Ox, Oy, hx, hy,          // scalars
           nx, ny, NP,
           Xp, Yp,                  // input
           dt, vx, vy, 
           F1, F2, F3, F4,
           Fup1, Fup2, Fup3, Fup4   // output
           );
}