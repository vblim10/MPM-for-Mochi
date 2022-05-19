/*==========================================================
 * DV_comp_2d.c - grid to particle interpolation
 *
 * Interpolates Eul. velocity v and takes its Jacobian (derivative)
 *        to return Lag. term DV
 *
 * The calling syntax is:
 *
 *		[DV1,DV2,DV3,DV4] = DV_comp_2d(Ox,Oy, hx,hy, nx,ny, Xp,Yp, vx,vy)
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
void DV_comp(const double Ox, const double Oy,                   // scalars
             const double hx, const double hy, 
             size_t nx, size_t ny, size_t NP,
             double *Xp, double *Yp,                             // input
             double *vx, double *vy,
             double *DV1, double *DV2, double *DV3, double *DV4  // output       
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
        DV1[pidx] = fp_curr1;
        DV2[pidx] = fp_curr2;
        DV3[pidx] = fp_curr3;
        DV4[pidx] = fp_curr4;
    }
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    // Declare variables
    double *Xp, *Yp;                    // Particle positions (NP x 1)
    double *vx, *vy;                    // grid vel v=(u+w)   (NG x 1)
    double *DV1, *DV2, *DV3, *DV4;      // DV                 (NP x 1)
    //const double Ox; const double Oy;   // "origin" of Eul. gridS
    //const double hx; const double hy;   // grid spacing
    size_t nx; size_t ny;               // Eul. grid partitions
    size_t NP;                          // # particles

    /* check for proper number of arguments */
    if (nrhs!=10) {
        mexErrMsgIdAndTxt("DV_comp_2d:nrhs","10 inputs required.");
    }
    if (nlhs!=4) {
        mexErrMsgIdAndTxt("DV_comp_2d:nlhs","4 outputs required.");
    }

    /* make sure the 6th input argument is type double */
    if (!mxIsDouble(prhs[6]) ||  mxIsComplex(prhs[6])) {
        mexErrMsgIdAndTxt("DV_comp_2d:notDouble","Xp must be type double.");
    }

   /* check that number of cols in input arguments is 1 */
    if (mxGetN(prhs[6])!=1 || mxGetN(prhs[7])!=1 ) {
        mexErrMsgIdAndTxt("DV_comp_2d:notRowVector","Xp, Yp must be a column vectors.");
    }
   
    /* check that number of rows in input arguments is 1 */
    if (mxGetM(prhs[8])!=1 || mxGetM(prhs[9])!=1  ) {
        mexErrMsgIdAndTxt("DV_comp_2d:notRowVector","vx, vy must be a row vectors.");
    }
   
    /* make sure the following input arguments are scalar */
    if (mxGetNumberOfElements(prhs[0])!=1 || mxGetNumberOfElements(prhs[1])!=1 ||
        mxGetNumberOfElements(prhs[2])!=1 || mxGetNumberOfElements(prhs[3])!=1 ||
        mxGetNumberOfElements(prhs[4])!=1 || mxGetNumberOfElements(prhs[5])!=1 ) {
        mexErrMsgIdAndTxt("DV_comp_2d:notScalar","Inputs Ox, Oy, hx, hy, nx, ny must be scalars.");
    }
   
    /* get the values of the scalar inputs  */
    const double Ox = mxGetScalar(prhs[0]);
    const double Oy = mxGetScalar(prhs[1]);
    const double hx = mxGetScalar(prhs[2]);
    const double hy = mxGetScalar(prhs[3]);
    nx = mxGetScalar(prhs[4]);
    ny = mxGetScalar(prhs[5]);

    /* create a pointer to the real data in Xp, Yp, vx, vy */
    #if MX_HAS_INTERLEAVED_COMPLEX      // get a C error if we don't do this if/else statement
    Xp = mxGetDoubles(prhs[6]);
    Yp = mxGetDoubles(prhs[7]);
    vx = mxGetDoubles(prhs[8]);
    vy = mxGetDoubles(prhs[9]);
    #else
    Xp = mxGetPr(prhs[6]);              // matlab recommends not to use mxGetPr
    Yp = mxGetPr(prhs[7]);
    vx = mxGetPr(prhs[8]);
    vy = mxGetPr(prhs[9]);
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
    DV1 = mxGetDoubles(plhs[0]);
    DV2 = mxGetDoubles(plhs[1]);
    DV3 = mxGetDoubles(plhs[2]);
    DV4 = mxGetDoubles(plhs[3]);
    #else
    DV1 = mxGetPr(plhs[0]); 
    DV2 = mxGetPr(plhs[1]);
    DV3 = mxGetPr(plhs[2]); 
    DV4 = mxGetPr(plhs[3]);
    #endif
   
    /* call the computational routine */
    DV_comp(Ox, Oy,                 // scalars
            hx, hy, 
            nx, ny, NP,
            Xp, Yp,                 // input
            vx, vy,
            DV1, DV2, DV3, DV4      // output       
            );
}