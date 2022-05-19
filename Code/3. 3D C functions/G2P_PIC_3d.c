/*==========================================================
 * G2P_PIC_3d.c - grid to particle interpolation
 *
 * Interpolates grid values fg
 *       to particle values fp
 *
 * The calling syntax is:
 *
 *		fp = G2P_PIC_3d(Ox,Oy,Oz, hx,hy,hz, nx,ny,nz, Xp,Yp,Zp, fg, bx,by,bz)
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
void G2P(const double Ox, const double Oy, const double Oz, // constants
         const double hx, const double hy, const double hz,
         size_t nx, size_t ny, size_t nz, size_t NP,
         double *Xp, double *Yp, double *Zp, double *fg,    // input
         double bx, double by, double bz,
         double *fp                                         // output       
        )
{
    // Constants
    const int d = 3;
   
    // Temp vars
    double fp_curr; 
    double xp_curr, yp_curr, zp_curr;
    double xg_curr, yg_curr, zg_curr, dx, dy, dz;
    int I, J, K;
    int i, j, k, gidx;
    double wh = pow(hx,bx)*pow(hy,by)*pow(hz,bz);   // need to divide by hx, hy, or hz 
                                                    // if doing grad. interp.
   
    // EVERYTHING CONDENSED IN THIS LOOP OVER PARTICLES
    for (int pidx=0; pidx<NP; pidx++) {
        // Reset temp vars
        fp_curr = 0.0;
        
        // Position for particle p
        xp_curr = Xp[pidx];
        yp_curr = Yp[pidx];
        zp_curr = Zp[pidx];
      
        // Index of reference grid node for particle p 
        I = floor( (xp_curr-Ox) / hx );  
        J = floor( (yp_curr-Oy) / hy );  
        K = floor( (zp_curr-Oz) / hz );
        
        // Loop over each level k
        for (int level=0; level<7; level++) {
            k = K + level - 3;          // k: Kp + (-3 to 3)
            // Check if level k is valid
            if(k>=0 && k<nz+1) {
                // Loop over stencil nodes
                for (int S=0; S<49; S++) {
                    i = floor(S/7);     // i: 0 to 6
                    j = S - (7*i);      // j: 0 to 6
                    i = I + i - 3;      // i: Ip + (-3 to 3)
                    j = J + j - 3;      // j: Jp + (-3 to 3)
                 
                    // Check if stencil indices are valid 
                    if (i>=0 && i<nx+1 && j>=0 && j<ny+1) {
                        // Compute gidx of (i,j,k) 
                        gidx = i*(ny+1) + j + k*(nx+1)*(ny+1); 
                        xg_curr = Ox + hx*i;
                        yg_curr = Oy + hy*j;
                        zg_curr = Oz + hz*k;
                    
                        // Compute Interpolation Weight
                        dx = (xp_curr - xg_curr)/hx; 
                        dy = (yp_curr - yg_curr)/hy; 
                        dz = (zp_curr - zg_curr)/hz;
                    
                        // PIC transfer
                        fp_curr += fg[gidx]*( N_weight2(dx,bx) * N_weight2(dy,by) * N_weight2(dz,bz) )/ wh; 
                    }
                }
            }
        }
        fp[pidx] = fp_curr;
    }
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    // Declare variables
    double *Xp; double *Yp; double *Zp; // Particle positions (NP x 1)
    double *fg;                         // grid values        (1 x NG)
    double *fp;                         // particle values    (1 x NP)
    double bx, by, bz;                  // indicator values (mult. by N, Nx, Ny, or Nz)
                                        // (bx,by,by) = (0,0,0), then N
                                        //            = (1,0,0), then Nx
                                        //            = (0,1,0), then Ny
                                        //            = (0,0,1), then Nz
    //const double Ox; const double Oy; const double Oz;  // "origin" of Eul. grid
    //const double hx; const double hy; const double hz;  // grid spacing
    size_t nx; size_t ny; size_t nz;    // Eul. grid partitions
    size_t NP;                          // # particles

    /* check for proper number of arguments */
    if (nrhs!=16) {
        mexErrMsgIdAndTxt("G2P_PIC_3d:nrhs","16 inputs required.");
    }
    if (nlhs!=1) {
        mexErrMsgIdAndTxt("G2P_PIC_3d:nlhs","One output required.");
    }

    /* make sure the 10th input argument is type double */
    if (!mxIsDouble(prhs[9]) ||  mxIsComplex(prhs[9])) {
        mexErrMsgIdAndTxt("G2P_PIC_3d:notDouble","Xp must be type double.");
    }

    /* check that number of cols in input arguments is 1 */
    if (mxGetN(prhs[9])!=1 || mxGetN(prhs[10])!=1 || mxGetN(prhs[11])!=1) {
        mexErrMsgIdAndTxt("G2P_PIC_3d:notColVector","Xp, Yp, Zp must be column vectors.");
    }
   
    /* check that number of rows in input arguments is 1 */
    if (mxGetM(prhs[12])!=1) {
        mexErrMsgIdAndTxt("G2P_PIC_3d:notRowVector","fg must be a row vector.");
    }
   
    /* make sure the following input arguments are scalar */
    if (mxGetNumberOfElements(prhs[0])!=1  || mxGetNumberOfElements(prhs[1])!=1  ||
        mxGetNumberOfElements(prhs[2])!=1  || mxGetNumberOfElements(prhs[3])!=1  ||
        mxGetNumberOfElements(prhs[4])!=1  || mxGetNumberOfElements(prhs[5])!=1  ||
        mxGetNumberOfElements(prhs[6])!=1  || mxGetNumberOfElements(prhs[7])!=1  ||
        mxGetNumberOfElements(prhs[8])!=1  || mxGetNumberOfElements(prhs[13])!=1 ||
        mxGetNumberOfElements(prhs[14])!=1 || mxGetNumberOfElements(prhs[15])!=1 ) {
        mexErrMsgIdAndTxt("G2P_PIC_3d:notScalar","Inputs Ox, Oy, Oz, hx, hy, hz, nx, ny, nz, bx, by, bz must be scalars.");
    }
   
    /* get the values of the scalar inputs  */
    const double Ox = mxGetScalar(prhs[0]);
    const double Oy = mxGetScalar(prhs[1]);
    const double Oz = mxGetScalar(prhs[2]);
    const double hx = mxGetScalar(prhs[3]);
    const double hy = mxGetScalar(prhs[4]);
    const double hz = mxGetScalar(prhs[5]);
    nx = mxGetScalar(prhs[6]);
    ny = mxGetScalar(prhs[7]);
    nz = mxGetScalar(prhs[8]);
    bx = mxGetScalar(prhs[13]);
    by = mxGetScalar(prhs[14]);
    bz = mxGetScalar(prhs[15]);
   
    /* create a pointer to the real data in Xp, Yp, Zp, fg */
    #if MX_HAS_INTERLEAVED_COMPLEX      // get a C error if we don't do this if/else statement
    Xp = mxGetDoubles(prhs[9]);
    Yp = mxGetDoubles(prhs[10]);
    Zp = mxGetDoubles(prhs[11]);
    fg = mxGetDoubles(prhs[12]);
    #else
    Xp = mxGetPr(prhs[9]);              // matlab recommends not to use mxGetPr
    Yp = mxGetPr(prhs[10]);
    Zp = mxGetPr(prhs[11]);
    fg = mxGetPr(prhs[12]);
    #endif
   
    /* get dimensions of the input matrix */
    NP = mxGetM(prhs[9]);

    /* create the output matrix */
    plhs[0] = mxCreateDoubleMatrix(1,(mwSize)NP,mxREAL);

    /* get a pointer to the real data in the output matrix */
    #if MX_HAS_INTERLEAVED_COMPLEX
    fp = mxGetDoubles(plhs[0]);
    #else
    fp = mxGetPr(plhs[0]);
    #endif
   
    /* call the computational routine */
    G2P(Ox, Oy, Oz,     // constants
        hx, hy, hz,
        nx, ny, nz, NP,
        Xp, Yp, Zp, fg, // input
        bx, by, bz,
        fp              // output       
        );
}