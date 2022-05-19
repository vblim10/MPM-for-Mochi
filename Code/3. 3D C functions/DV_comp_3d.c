/*==========================================================
 * DV_comp_3d.c - grid to particle interpolation
 *
 * Interpolates Eul. velocity v and takes its Jacobian (derivative)
 *        to return Lag. term DV
 *
 * The calling syntax is:
 *
 *		[DV1, DV2, DV3, DV4, DV5, DV6, DV7, DV8, DV9] 
 *      = DV_comp_3d(Ox,Oy,Oz, hx,hy,hz, nx,ny,nz, Xp,Yp,Zp, vx,vy,vz)
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

// Deformation gradient computation routine
void DV_comp(const double Ox, const double Oy, const double Oz, // constants
             const double hx, const double hy, const double hz,
             size_t nx, size_t ny, size_t nz, size_t NP,
             double *Xp, double *Yp, double *Zp,                // input
             double *vx, double *vy, double *vz,
             double *DV1, double *DV2, double *DV3,             // output 
             double *DV4, double *DV5, double *DV6,
             double *DV7, double *DV8, double *DV9
             )
{
    // Constants
    const int d = 3;
   
    // Temp vars
    double fp1_curr, fp2_curr, fp3_curr; 
    double fp4_curr, fp5_curr, fp6_curr; 
    double fp7_curr, fp8_curr, fp9_curr; 
    double xp_curr, yp_curr, zp_curr;
    double xg_curr, yg_curr, zg_curr, dx, dy, dz;
    int I, J, K;
    int i, j, k, gidx;
    double Nx, Ny, Nz;
   
    // EVERYTHING CONDENSED IN THIS LOOP OVER PARTICLES
    for (int pidx=0; pidx<NP; pidx++) {
        // Reset temp vars
        fp1_curr = 0.0; fp4_curr = 0.0; fp7_curr = 0.0;
        fp2_curr = 0.0; fp5_curr = 0.0; fp8_curr = 0.0;
        fp3_curr = 0.0; fp6_curr = 0.0; fp9_curr = 0.0;
        
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
                        Nx = ( N_weight2(dx,1) * N_weight2(dy,0) * N_weight2(dz,0) ) / hx;
                        Ny = ( N_weight2(dx,0) * N_weight2(dy,1) * N_weight2(dz,0) ) / hy;
                        Nz = ( N_weight2(dx,0) * N_weight2(dy,0) * N_weight2(dz,1) ) / hz;
                    
                        // PIC transfer
                        fp1_curr += vx[gidx] * Nx; 
                        fp2_curr += vy[gidx] * Nx;
                        fp3_curr += vz[gidx] * Nx;
                        
                        fp4_curr += vx[gidx] * Ny; 
                        fp5_curr += vy[gidx] * Ny;
                        fp6_curr += vz[gidx] * Ny;
                        
                        fp7_curr += vx[gidx] * Nz; 
                        fp8_curr += vy[gidx] * Nz;
                        fp9_curr += vz[gidx] * Nz;
                    }
                }
            }
        }        
        // DV
        DV1[pidx] = fp1_curr;
        DV2[pidx] = fp2_curr;
        DV3[pidx] = fp3_curr;
        DV4[pidx] = fp4_curr;
        DV5[pidx] = fp5_curr;
        DV6[pidx] = fp6_curr;
        DV7[pidx] = fp7_curr;
        DV8[pidx] = fp8_curr;
        DV9[pidx] = fp9_curr;
    }
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    // Declare variables
    double *Xp; double *Yp; double *Zp; // Particle positions (NP x 1)
    double *vx,   *vy,   *vz;           // grid vel v=(u+w)   (NG x 1)
    double *DV1, *DV4, *DV7;                   
    double *DV2, *DV5, *DV8;            // DV                 (NP x 1)
    double *DV3, *DV6, *DV9;
    //const double Ox; const double Oy; const double Oz;  // "origin" of Eul. grid
    //const double hx; const double hy; const double hz;  // grid spacing
    size_t nx; size_t ny; size_t nz;    // Eul. grid partitions
    size_t NP;                          // # particles

    /* check for proper number of arguments */
    if (nrhs!=15) {
        mexErrMsgIdAndTxt("DV_comp_3d:nrhs","15 inputs required.");
    }
    if (nlhs!=9) {
        mexErrMsgIdAndTxt("DV_comp_3d:nlhs","9 outputs required.");
    }

    /* make sure the 10th input argument is type double */
    if ( !mxIsDouble(prhs[9]) ||  mxIsComplex(prhs[9])) {
        mexErrMsgIdAndTxt("DV_comp_3d:notDouble","Xp must be type double.");
    }

    /* check that number of cols in input arguments is 1 */
    if (mxGetN(prhs[9])!=1 || mxGetN(prhs[10])!=1 || mxGetN(prhs[11])!=1) {
        mexErrMsgIdAndTxt("DV_comp_3d:notColVector","Xp, Yp, Zp must be column vectors.");
    }
   
    /* check that number of rows in input arguments is 1 */
    if (mxGetM(prhs[12])!=1 || mxGetM(prhs[13])!=1 || mxGetM(prhs[14])!=1 ) {
        mexErrMsgIdAndTxt("DV_comp_3d:notRowVector","vx, vy, vz must be row vectors.");
    }
   
    /* make sure the following input arguments are scalar */
    if (mxGetNumberOfElements(prhs[0])!=1  || mxGetNumberOfElements(prhs[1])!=1  ||
        mxGetNumberOfElements(prhs[2])!=1  || mxGetNumberOfElements(prhs[3])!=1  ||
        mxGetNumberOfElements(prhs[4])!=1  || mxGetNumberOfElements(prhs[5])!=1  ||
        mxGetNumberOfElements(prhs[6])!=1  || mxGetNumberOfElements(prhs[7])!=1  ||
        mxGetNumberOfElements(prhs[8])!=1  ) {
        mexErrMsgIdAndTxt("DV_comp_3d:notScalar","Inputs Ox, Oy, Oz, hx, hy, hz, nx, ny, nz must be scalars.");
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
   
    /* create a pointer to the real data in Xp, Yp, Zp, vx, vy, vz */
    #if MX_HAS_INTERLEAVED_COMPLEX      // get a C error if we don't do this if/else statement
    Xp = mxGetDoubles(prhs[9]);
    Yp = mxGetDoubles(prhs[10]);
    Zp = mxGetDoubles(prhs[11]);
    vx = mxGetDoubles(prhs[12]);
    vy = mxGetDoubles(prhs[13]);
    vz = mxGetDoubles(prhs[14]);
    #else
    Xp = mxGetPr(prhs[9]);
    Yp = mxGetPr(prhs[10]);
    Zp = mxGetPr(prhs[11]);
    vx = mxGetPr(prhs[12]);
    vy = mxGetPr(prhs[13]);
    vz = mxGetPr(prhs[14]);
    #endif
   
    /* get dimensions of the input matrix */
    NP = mxGetM(prhs[9]);

    /* create the output matrices */
    plhs[0] = mxCreateDoubleMatrix(1,(mwSize)NP,mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1,(mwSize)NP,mxREAL);
    plhs[2] = mxCreateDoubleMatrix(1,(mwSize)NP,mxREAL);
    plhs[3] = mxCreateDoubleMatrix(1,(mwSize)NP,mxREAL);
    plhs[4] = mxCreateDoubleMatrix(1,(mwSize)NP,mxREAL);
    plhs[5] = mxCreateDoubleMatrix(1,(mwSize)NP,mxREAL);
    plhs[6] = mxCreateDoubleMatrix(1,(mwSize)NP,mxREAL);
    plhs[7] = mxCreateDoubleMatrix(1,(mwSize)NP,mxREAL);
    plhs[8] = mxCreateDoubleMatrix(1,(mwSize)NP,mxREAL);

    /* get a pointer to the real data in the output matrix */
    #if MX_HAS_INTERLEAVED_COMPLEX
    DV1 = mxGetDoubles(plhs[0]);
    DV2 = mxGetDoubles(plhs[1]);
    DV3 = mxGetDoubles(plhs[2]);
    DV4 = mxGetDoubles(plhs[3]);
    DV5 = mxGetDoubles(plhs[4]);
    DV6 = mxGetDoubles(plhs[5]);
    DV7 = mxGetDoubles(plhs[6]);
    DV8 = mxGetDoubles(plhs[7]);
    DV9 = mxGetDoubles(plhs[8]);
    #else
    DV1 = mxGetPr(plhs[0]);
    DV2 = mxGetPr(plhs[1]);
    DV3 = mxGetPr(plhs[2]);
    DV4 = mxGetPr(plhs[3]);
    DV5 = mxGetPr(plhs[4]);
    DV6 = mxGetPr(plhs[5]);
    DV7 = mxGetPr(plhs[6]);
    DV8 = mxGetPr(plhs[7]);
    DV9 = mxGetPr(plhs[8]);
    #endif
   
    /* call the computational routine */
    DV_comp(Ox, Oy, Oz,     // constants
            hx, hy, hz,
            nx, ny, nz, NP,
            Xp, Yp, Zp,     // input
            vx, vy, vz,
            DV1, DV2, DV3,  // output 
            DV4, DV5, DV6,
            DV7, DV8, DV9
            );
}