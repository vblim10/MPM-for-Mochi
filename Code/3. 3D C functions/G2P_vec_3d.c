/*==========================================================
 * G2P_vec_3d.c - grid to particle interpolation
 *
 * Interpolates grid values fg1, fg2, fg3
 *       to particle values fp1, fp2, fp3
 *
 * The calling syntax is:
 *
 *		[fp1,fp2,fp3] = G2P_vec_3d(Ox,Oy,Oz, hx,hy,hz, nx,ny,nz, Xp,Yp,Zp, fg1,fg2,fg3, bx,by,bz)
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
         double *Xp, double *Yp, double *Zp,                // input
         double *fg1, double *fg2, double *fg3,
         double bx, double by, double bz,
         double *fp1, double *fp2, double *fp3              // output       
        )
{
    // Constants
    const int d = 3;
   
    // Temp vars
    double fp1_curr, fp2_curr, fp3_curr; 
    double xp_curr, yp_curr, zp_curr;
    double xg_curr, yg_curr, zg_curr, dx, dy, dz;
    int I, J, K;
    int i, j, k, gidx;
    double w = 0.0;
    double wh = pow(hx,bx)*pow(hy,by)*pow(hz,bz);   // need to divide by hx, hy, or hz 
                                                    // if doing grad. interp.
   
    // EVERYTHING CONDENSED IN THIS LOOP OVER PARTICLES
    for (int pidx=0; pidx<NP; pidx++) {
        // Reset temp vars
        fp1_curr = 0.0;
        fp2_curr = 0.0; 
        fp3_curr = 0.0;
        
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
                        w = ( N_weight2(dx,bx) * N_weight2(dy,by) * N_weight2(dz,bz) )/ wh; 
                    
                        // PIC transfer
                        fp1_curr += fg1[gidx] * w; 
                        fp2_curr += fg2[gidx] * w;
                        fp3_curr += fg3[gidx] * w;
                    }
                }
            }
        }
        fp1[pidx] = fp1_curr;
        fp2[pidx] = fp2_curr;
        fp3[pidx] = fp3_curr;
    }
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    // Declare variables
    double *Xp; double *Yp; double *Zp;     // Particle positions (NP x 1)
    double *fg1; double *fg2; double *fg3;  // grid values        (1 x NG)
    double *fp1; double *fp2; double *fp3;  // particle values    (1 x NP)
    double bx, by, bz;  // indicator values (mult. by N, Nx, Ny, or Nz)
                        // (bx,by,by) = (0,0,0), then N
                        //            = (1,0,0), then Nx
                        //            = (0,1,0), then Ny
                        //            = (0,0,1), then Nz
    //const double Ox; const double Oy; const double Oz;  // "origin" of Eul. grid
    //const double hx; const double hy; const double hz;  // grid spacing
    size_t nx; size_t ny; size_t nz;    // Eul. grid partitions
    size_t NP;                          // # particles

    /* check for proper number of arguments */
    if (nrhs!=18) {
        mexErrMsgIdAndTxt("G2P_vec_3d:nrhs","18 inputs required.");
    }
    if (nlhs!=3) {
        mexErrMsgIdAndTxt("G2P_vec_3d:nlhs","3 outputs required.");
    }

    /* make sure the 10th input argument is type double */
    if (!mxIsDouble(prhs[9]) ||  mxIsComplex(prhs[9])) {
        mexErrMsgIdAndTxt("G2P_vec_3d:notDouble","Xp must be type double.");
    }

    /* check that number of cols in input arguments is 1 */
    if (mxGetN(prhs[9])!=1 || mxGetN(prhs[10])!=1 || mxGetN(prhs[11])!=1) {
        mexErrMsgIdAndTxt("G2P_vec_3d:notColVector","Xp, Yp, Zp must be column vectors.");
    }
   
    /* check that number of rows in input arguments is 1 */
    if (mxGetM(prhs[12])!=1 || mxGetM(prhs[13])!=1 || mxGetM(prhs[14])!=1) {
        mexErrMsgIdAndTxt("G2P_vec_3d:notRowVector","fg1, fg2, fg3 must be row vectors.");
    }
   
    /* make sure the following input arguments are scalar */
    if (mxGetNumberOfElements(prhs[0])!=1  || mxGetNumberOfElements(prhs[1])!=1  ||
        mxGetNumberOfElements(prhs[2])!=1  || mxGetNumberOfElements(prhs[3])!=1  ||
        mxGetNumberOfElements(prhs[4])!=1  || mxGetNumberOfElements(prhs[5])!=1  ||
        mxGetNumberOfElements(prhs[6])!=1  || mxGetNumberOfElements(prhs[7])!=1  ||
        mxGetNumberOfElements(prhs[8])!=1  || mxGetNumberOfElements(prhs[15])!=1 ||
        mxGetNumberOfElements(prhs[16])!=1 || mxGetNumberOfElements(prhs[17])!=1 ) {
        mexErrMsgIdAndTxt("G2P_vec_3d:notScalar","Inputs Ox, Oy, Oz, hx, hy, hz, nx, ny, nz, bx, by, bz must be scalars.");
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
    bx = mxGetScalar(prhs[15]);
    by = mxGetScalar(prhs[16]);
    bz = mxGetScalar(prhs[17]);
   
    /* create a pointer to the real data in Xp, Yp, Zp, fg1, fg2, fg3 */
    #if MX_HAS_INTERLEAVED_COMPLEX      // get a C error if we don't do this if/else statement
    Xp = mxGetDoubles(prhs[9]);
    Yp = mxGetDoubles(prhs[10]);
    Zp = mxGetDoubles(prhs[11]);
    fg1 = mxGetDoubles(prhs[12]);
    fg2 = mxGetDoubles(prhs[13]);
    fg3 = mxGetDoubles(prhs[14]);
    #else
    Xp = mxGetPr(prhs[9]);              // matlab recommends not to use mxGetPr
    Yp = mxGetPr(prhs[10]);
    Zp = mxGetPr(prhs[11]);
    fg1 = mxGetPr(prhs[12]);
    fg2 = mxGetPr(prhs[13]);
    fg3 = mxGetPr(prhs[14]);
    #endif
   
    /* get dimensions of the input matrix */
    NP = mxGetM(prhs[9]);

    /* create the output matrices */
    plhs[0] = mxCreateDoubleMatrix(1,(mwSize)NP,mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1,(mwSize)NP,mxREAL);
    plhs[2] = mxCreateDoubleMatrix(1,(mwSize)NP,mxREAL);

    /* get a pointer to the real data in the output matrix */
    #if MX_HAS_INTERLEAVED_COMPLEX
    fp1 = mxGetDoubles(plhs[0]);
    fp2 = mxGetDoubles(plhs[1]);
    fp3 = mxGetDoubles(plhs[2]);
    #else
    fp1 = mxGetPr(plhs[0]);
    fp2 = mxGetPr(plhs[1]);
    fp3 = mxGetPr(plhs[2]);
    #endif
   
    /* call the computational routine */
    G2P(Ox, Oy, Oz,     // constants
        hx, hy, hz,
        nx, ny, nz, NP,
        Xp, Yp, Zp,     // input
        fg1, fg2, fg3,
        bx, by, bz,
        fp1, fp2, fp3   // output       
        );
}