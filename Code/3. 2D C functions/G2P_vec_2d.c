/*==========================================================
 * G2P_vec_2d.c - grid to particle interpolation
 *
 * Interpolates grid values fg1, fg2
 *       to particle values fp1, fp2
 *
 * The calling syntax is:
 *
 *		[fp1,fp2] = G2P_vec_2d(Ox,Oy, hx,hy, nx,ny, Xp,Yp, fg1,fg2, bx,by)
 *
 * This is a MEX-file for MATLAB.
 *
 *========================================================*/

#include "mex.h"
#include <stdio.h>
#include <math.h>

// Weight Function
double N_weight(double x) // old version
{
    x = fabs(x);    // *** abs function is for ints only (in C) ***
    double N = 0.0;
    if (x<0.5)
        N += 0.75 - x*x;
    else if (x<1.5)
        N += 0.5*(1.5-x)*(1.5-x);

    return N;
}

double N_weight2(double x, double bw) // new version
{
    // bw = 1: returns Nx 
    // bw = 0: returns N 
    int sgn = (x>0)?1:-1;
   
    x = fabs(x);    // *** abs function is for ints only (in C) ***
    double N = 0.0;
    if (x<0.5)
        N += (1-bw)*(0.75 - x*x)        + bw*sgn*(-2*x);
    else if (x<1.5)
        N += (1-bw)*0.5*(1.5-x)*(1.5-x) + bw*sgn*(-1.5+x);

    return N;
}

// G2P interpolation routine
void G2P(const double Ox, const double Oy,      // constants
         const double hx, const double hy, 
         size_t nx, size_t ny, size_t NP,
         double *Xp, double *Yp,                // input
         double *fg1, double *fg2,
         double bx, double by,
         double *fp1, double *fp2               // output       
         )
{
   // Constants
   const int d = 2;
   
   // Temp vars
   double fp1_curr, fp2_curr; 
   double xp_curr, yp_curr;
   double xg_curr, yg_curr, dx, dy;
   int I, J;
   int l, k, gidx;
   double w = 0.0;
   double wh = pow(hx,bx)*pow(hy,by);  // need to divide by hx or hy if doing grad. interp.
   
   
    // EVERYTHING CONDENSED IN THIS LOOP OVER PARTICLES
    for (int pidx=0; pidx<NP; pidx++) {
        fp1_curr = 0.0;
        fp2_curr = 0.0;
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
                w = ( N_weight2(dx,bx) * N_weight2(dy,by) )/ wh;
            
                // PIC transfer
                fp1_curr += fg1[gidx] * w;
                fp2_curr += fg2[gidx] * w;
            }
        }
        fp1[pidx] = fp1_curr;
        fp2[pidx] = fp2_curr;
    }
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    // Declare variables
    double *Xp, *Yp;        // Particle positions (NP x 1)
    double *fg1, *fg2;      // grid values        (1 x NG)
    double *fp1, *fp2;      // particle values    (1 x NP)
    double bx, by;          // indicator values (mult. by N, Nx, or Ny)
                            // (bx,by) = (0,0), then N
                            //         = (1,0), then Nx
                            //         = (0,1), then Ny
    //const double Ox; const double Oy;   // "origin" of Eul. grid
    //const double hx; const double hy;   // grid spacing
    size_t nx; size_t ny;   // Eul. grid partitions
    size_t NP;              // # particles

   /* check for proper number of arguments */
    if (nrhs!=12) {
        mexErrMsgIdAndTxt("G2P_vec_2d:nrhs","12 inputs required.");
    }
    if (nlhs!=2) {
        mexErrMsgIdAndTxt("G2P_vec_2d:nlhs","2 outputs required.");
    }

    /* make sure the 7th input argument is type double */
    if (!mxIsDouble(prhs[6]) ||  mxIsComplex(prhs[6])) {
        mexErrMsgIdAndTxt("G2P_vec_2d:notDouble","Xp must be type double.");
    }

    /* check that number of cols in input arguments is 1 */
    if (mxGetN(prhs[6])!=1 || mxGetN(prhs[7])!=1 ) {
        mexErrMsgIdAndTxt("G2P_vec_2d:notColVector","Xp, Yp must be a column vectors.");
    }
   
    /* check that number of rows in input arguments is 1 */
    if (mxGetM(prhs[8])!=1 || mxGetM(prhs[9])!=1 ) {
        mexErrMsgIdAndTxt("G2P_vec_2d:notRowVector","fg1, fg2 must be a row vectors.");
    }
   
    /* make sure the following input arguments are scalar */
    if (mxGetNumberOfElements(prhs[0])!=1 || mxGetNumberOfElements(prhs[1])!=1 ||
        mxGetNumberOfElements(prhs[2])!=1 || mxGetNumberOfElements(prhs[3])!=1 ||
        mxGetNumberOfElements(prhs[4])!=1 || mxGetNumberOfElements(prhs[5])!=1 ||
        mxGetNumberOfElements(prhs[10])!=1 || mxGetNumberOfElements(prhs[11])!=1 ) {
        mexErrMsgIdAndTxt("G2P_vec_2d:notScalar","Inputs Ox, Oy, hx, hy, nx, ny, bx, by must be scalars.");
    }

    /* get the values of the scalar inputs  */
    const double Ox = mxGetScalar(prhs[0]);
    const double Oy = mxGetScalar(prhs[1]);
    const double hx = mxGetScalar(prhs[2]);
    const double hy = mxGetScalar(prhs[3]);
    nx = mxGetScalar(prhs[4]);
    ny = mxGetScalar(prhs[5]);
    bx = mxGetScalar(prhs[10]);
    by = mxGetScalar(prhs[11]);
   
    /* create a pointer to the real data in Xp, Yp, fg1, fg2 */
    #if MX_HAS_INTERLEAVED_COMPLEX      // get a C error if we don't do this if/else statement
    Xp = mxGetDoubles(prhs[6]);
    Yp = mxGetDoubles(prhs[7]);
    fg1 = mxGetDoubles(prhs[8]);
    fg2 = mxGetDoubles(prhs[9]);
    #else
    Xp = mxGetPr(prhs[6]);              // matlab recommends not to use mxGetPr
    Yp = mxGetPr(prhs[7]);
    fg1 = mxGetPr(prhs[8]);
    fg2 = mxGetPr(prhs[9]);
    #endif
   
    /* get dimensions of the input matrix */
    NP = mxGetM(prhs[6]);

    /* create the output matrices */
    plhs[0] = mxCreateDoubleMatrix(1,(mwSize)NP,mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1,(mwSize)NP,mxREAL);

    /* get a pointer to the real data in the output matrix */
    #if MX_HAS_INTERLEAVED_COMPLEX
    fp1 = mxGetDoubles(plhs[0]);
    fp2 = mxGetDoubles(plhs[1]);
    #else
    fp1 = mxGetPr(plhs[0]);
    fp2 = mxGetPr(plhs[1]);
    #endif
   
    /* call the computational routine */
    G2P(Ox, Oy, hx, hy, nx, ny, NP,  // scalars
        Xp, Yp,                      // input
        fg1, fg2, bx, by,
        fp1, fp2                     // output
        );
}