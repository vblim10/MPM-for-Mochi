/*==========================================================
 * G2P_APIC_3d.c - grid to particle interpolation
 *
 * 1. Interpolates grid values (xup1,xup2,xup3) and (w1,w2,w3)
 *          to particle values (Phi1,Phi2,Phi3) and (Vp1,Vp2,Vp3)
 *
 * 2. Updates APIC var Bp
 *
 * The calling syntax is:
 *
 *		[Phi1,Phi2,Phi3, Vp1,Vp2,Vp3, B1,B2,B3,B4,B5,B6,B7,B8,B9] 
 *      = G2P_APIC_3d(Ox,Oy,Oz, hx,hy,hz, nx,ny,nz, Xp,Yp,Zp, xup1,xup2,xup3, w1,w2,w3)
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
         double *xup1, double *xup2, double *xup3, 
         double *w1, double *w2, double *w3,  
         double *Phi1, double *Phi2, double *Phi3,          // output  
         double *Vp1, double *Vp2, double *Vp3,
         double *B1, double *B2, double *B3, 
         double *B4, double *B5, double *B6, 
         double *B7, double *B8, double *B9      
        )
{
    // Constants
    const int d = 3;
   
    // Temp vars
    double fp1_curr, fp2_curr, fp3_curr;    // for Phi update
    double fp4_curr, fp5_curr, fp6_curr;    // for Vp update
    double xp_curr, yp_curr, zp_curr;
    double xg_curr, yg_curr, zg_curr, dx, dy, dz;
    int I, J, K;
    int i, j, k, gidx;
    //double N[7][49];
    double bl1, bl2, bl3; 
    double br1, br2, br3;
   
    // EVERYTHING CONDENSED IN THIS LOOP OVER PARTICLES
    for (int pidx=0; pidx<NP; pidx++) {
        // Reset temp vars
        double N[7][49] = {0.0};    // array of all zeros 
        fp1_curr = 0.0; fp2_curr = 0.0; fp3_curr = 0.0;
        fp4_curr = 0.0; fp5_curr = 0.0; fp6_curr = 0.0;
        
        // Position for particle p
        xp_curr = Xp[pidx];
        yp_curr = Yp[pidx];
        zp_curr = Zp[pidx];
      
        // Index of reference grid node for particle p 
        I = floor( (xp_curr-Ox) / hx );  
        J = floor( (yp_curr-Oy) / hy );  
        K = floor( (zp_curr-Oz) / hz );
        
        // 1st (*) loop
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
                        N[level][S] = ( N_weight2(dx,0) * N_weight2(dy,0) * N_weight2(dz,0) ); 
                        
                        // PIC transfer
                        fp1_curr += xup1[gidx] * N[level][S]; 
                        fp2_curr += xup2[gidx] * N[level][S];
                        fp3_curr += xup3[gidx] * N[level][S];
                        fp4_curr += w1[gidx] * N[level][S];
                        fp5_curr += w2[gidx] * N[level][S];
                        fp6_curr += w3[gidx] * N[level][S];
                    }
                }
            }
        }
        Phi1[pidx] = fp1_curr;
        Phi2[pidx] = fp2_curr;
        Phi3[pidx] = fp3_curr;
        Vp1[pidx] = fp4_curr;
        Vp2[pidx] = fp5_curr;
        Vp3[pidx] = fp6_curr;
        
        // 2nd (*) loop (need a 2nd double-loop b/c we need Phi)
        for (int level=0; level<7; level++) {
            k = K + level - 3;          // k: Kp + (-3 to 3)
            
            // Loop over stencil nodes
            for (int S=0; S<49; S++) {
                
                if (N[level][S] > 0.0) {// don't need to "validate" indices
                    i = floor(S/7);     // i: 0 to 6
                    j = S - (7*i);      // j: 0 to 6
                    i = I + i - 3;      // i: Ip + (-3 to 3)
                    j = J + j - 3;      // j: Jp + (-3 to 3)    
                    
                    // Compute gidx of (i,j,k) 
                    gidx = i*(ny+1) + j + k*(nx+1)*(ny+1); 
                    xg_curr = Ox + hx*i;
                    yg_curr = Oy + hy*j;
                    zg_curr = Oz + hz*k;
                    
                    bl1 = xg_curr - xp_curr + xup1[gidx] - Phi1[pidx];
                    bl2 = yg_curr - yp_curr + xup2[gidx] - Phi2[pidx];
                    bl3 = zg_curr - zp_curr + xup3[gidx] - Phi3[pidx];
                    
                    br1 = xg_curr - xp_curr - xup1[gidx] + Phi1[pidx];
                    br2 = yg_curr - yp_curr - xup2[gidx] + Phi2[pidx];
                    br3 = zg_curr - zp_curr - xup3[gidx] + Phi3[pidx];
                    
                    // Compute APIC var Bp
                    B1[pidx] += 0.5 * ( w1[gidx]*bl1 + br1*w1[gidx] ) * N[level][S];
                    B2[pidx] += 0.5 * ( w2[gidx]*bl1 + br2*w1[gidx] ) * N[level][S];
                    B3[pidx] += 0.5 * ( w3[gidx]*bl1 + br3*w1[gidx] ) * N[level][S];
                    
                    B4[pidx] += 0.5 * ( w1[gidx]*bl2 + br1*w2[gidx] ) * N[level][S];
                    B5[pidx] += 0.5 * ( w2[gidx]*bl2 + br2*w2[gidx] ) * N[level][S];
                    B6[pidx] += 0.5 * ( w3[gidx]*bl2 + br3*w2[gidx] ) * N[level][S];
                    
                    B7[pidx] += 0.5 * ( w1[gidx]*bl3 + br1*w3[gidx] ) * N[level][S];
                    B8[pidx] += 0.5 * ( w2[gidx]*bl3 + br2*w3[gidx] ) * N[level][S];
                    B9[pidx] += 0.5 * ( w3[gidx]*bl3 + br3*w3[gidx] ) * N[level][S];
                }
            }
        }
    }
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    // Declare variables
    double *Xp; double *Yp; double *Zp; // Particle positions (NP x 1)
    double *xup1, *xup2, *xup3;         // grid values        (1 x NG)
    double *w1, *w2, *w3;               //                    (1 x NG)
    double *Phi1, *Phi2, *Phi3;         // particle values    (NP x 1)
    double *Vp1, *Vp2, *Vp3;            //                    (1 x NP)
    double *B1, *B2, *B3;               // APIC var Bp        (1 x NP)
    double *B4, *B5, *B6;               //                    (1 x NP)
    double *B7, *B8, *B9;               //                    (1 x NP)
    //const double Ox; const double Oy; const double Oz;  // "origin" of Eul. grid
    //const double hx; const double hy; const double hz;  // grid spacing
    size_t nx; size_t ny; size_t nz;    // Eul. grid partitions
    size_t NP;                          // # particles

    /* check for proper number of arguments */
    if (nrhs!=18) {
        mexErrMsgIdAndTxt("G2P_APIC_3d:nrhs","18 inputs required.");
    }
    if (nlhs!=15) {
        mexErrMsgIdAndTxt("G2P_APIC_3d:nlhs","15 outputs required.");
    }

    /* make sure the 10th input argument is type double */
    if (!mxIsDouble(prhs[9]) ||  mxIsComplex(prhs[9])) {
        mexErrMsgIdAndTxt("G2P_APIC_3d:notDouble","Xp must be type double.");
    }

    /* check that number of cols in input arguments is 1 */
    if (mxGetN(prhs[9])!=1 || mxGetN(prhs[10])!=1 || mxGetN(prhs[11])!=1) {
        mexErrMsgIdAndTxt("G2P_APIC_3d:notColVector","Xp, Yp, Zp must be column vectors.");
    }
   
    /* check that number of rows in input arguments is 1 */
    if (mxGetM(prhs[12])!=1 || mxGetM(prhs[13])!=1 || mxGetM(prhs[14])!=1 || 
        mxGetM(prhs[15])!=1 || mxGetM(prhs[16])!=1 || mxGetM(prhs[17])!=1) {
        mexErrMsgIdAndTxt("G2P_APIC_3d:notRowVector","xup1, xup2, xup3, w1, w2, w3 must be row vectors.");
    }
   
    /* make sure the following input arguments are scalar */
    if (mxGetNumberOfElements(prhs[0])!=1  || mxGetNumberOfElements(prhs[1])!=1  ||
        mxGetNumberOfElements(prhs[2])!=1  || mxGetNumberOfElements(prhs[3])!=1  ||
        mxGetNumberOfElements(prhs[4])!=1  || mxGetNumberOfElements(prhs[5])!=1  ||
        mxGetNumberOfElements(prhs[6])!=1  || mxGetNumberOfElements(prhs[7])!=1  ||
        mxGetNumberOfElements(prhs[8])!=1  ) {
        mexErrMsgIdAndTxt("G2P_APIC_3d:notScalar","Inputs Ox, Oy, Oz, hx, hy, hz, nx, ny, nz must be scalars.");
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
   
    /* create a pointer to the real data in Xp, Yp, Zp, xup1, xup2, xup3, w1, w2, w3 */
    #if MX_HAS_INTERLEAVED_COMPLEX      // get a C error if we don't do this if/else statement
    Xp = mxGetDoubles(prhs[9]);
    Yp = mxGetDoubles(prhs[10]);
    Zp = mxGetDoubles(prhs[11]);
    xup1 = mxGetDoubles(prhs[12]);
    xup2 = mxGetDoubles(prhs[13]);
    xup3 = mxGetDoubles(prhs[14]);
    w1 = mxGetDoubles(prhs[15]);
    w2 = mxGetDoubles(prhs[16]);
    w3 = mxGetDoubles(prhs[17]);
    #else
    Xp = mxGetPr(prhs[9]);              // matlab recommends not to use mxGetPr
    Yp = mxGetPr(prhs[10]);
    Zp = mxGetPr(prhs[11]);
    xup1 = mxGetPr(prhs[12]);
    xup2 = mxGetPr(prhs[13]);
    xup3 = mxGetPr(prhs[14]);
    w1 = mxGetPr(prhs[15]);
    w2 = mxGetPr(prhs[16]);
    w3 = mxGetPr(prhs[17]);
    #endif
   
    /* get dimensions of the input matrix */
    NP = mxGetM(prhs[9]);

    /* create the output matrices */
    plhs[0] = mxCreateDoubleMatrix((mwSize)NP,1,mxREAL);
    plhs[1] = mxCreateDoubleMatrix((mwSize)NP,1,mxREAL);
    plhs[2] = mxCreateDoubleMatrix((mwSize)NP,1,mxREAL);
    plhs[3] = mxCreateDoubleMatrix(1,(mwSize)NP,mxREAL);
    plhs[4] = mxCreateDoubleMatrix(1,(mwSize)NP,mxREAL);
    plhs[5] = mxCreateDoubleMatrix(1,(mwSize)NP,mxREAL);
    plhs[6] = mxCreateDoubleMatrix(1,(mwSize)NP,mxREAL);
    plhs[7] = mxCreateDoubleMatrix(1,(mwSize)NP,mxREAL);
    plhs[8] = mxCreateDoubleMatrix(1,(mwSize)NP,mxREAL);
    plhs[9] = mxCreateDoubleMatrix(1,(mwSize)NP,mxREAL);
    plhs[10] = mxCreateDoubleMatrix(1,(mwSize)NP,mxREAL);
    plhs[11] = mxCreateDoubleMatrix(1,(mwSize)NP,mxREAL);
    plhs[12] = mxCreateDoubleMatrix(1,(mwSize)NP,mxREAL);
    plhs[13] = mxCreateDoubleMatrix(1,(mwSize)NP,mxREAL);
    plhs[14] = mxCreateDoubleMatrix(1,(mwSize)NP,mxREAL);

    /* get a pointer to the real data in the output matrix */
    #if MX_HAS_INTERLEAVED_COMPLEX
    Phi1 = mxGetDoubles(plhs[0]);
    Phi2 = mxGetDoubles(plhs[1]);
    Phi3 = mxGetDoubles(plhs[2]);
    Vp1 = mxGetDoubles(plhs[3]);
    Vp2 = mxGetDoubles(plhs[4]);
    Vp3 = mxGetDoubles(plhs[5]);
    B1 = mxGetDoubles(plhs[6]);
    B2 = mxGetDoubles(plhs[7]);
    B3 = mxGetDoubles(plhs[8]);
    B4 = mxGetDoubles(plhs[9]);
    B5 = mxGetDoubles(plhs[10]);
    B6 = mxGetDoubles(plhs[11]);
    B7 = mxGetDoubles(plhs[12]);
    B8 = mxGetDoubles(plhs[13]);
    B9 = mxGetDoubles(plhs[14]);
    #else
    Phi1 = mxGetPr(plhs[0]);
    Phi2 = mxGetPr(plhs[1]);
    Phi3 = mxGetPr(plhs[2]);
    Vp1 = mxGetPr(plhs[3]);
    Vp2 = mxGetPr(plhs[4]);
    Vp3 = mxGetPr(plhs[5]);
    B1 = mxGetPr(plhs[6]);
    B2 = mxGetPr(plhs[7]);
    B3 = mxGetPr(plhs[8]);
    B4 = mxGetPr(plhs[9]);
    B5 = mxGetPr(plhs[10]);
    B6 = mxGetPr(plhs[11]);
    B7 = mxGetPr(plhs[12]);
    B8 = mxGetPr(plhs[13]);
    B9 = mxGetPr(plhs[14]);
    #endif
   
    /* call the computational routine */
    G2P(Ox, Oy, Oz,         // constants
        hx, hy, hz,
        nx, ny, nz, NP,
        Xp, Yp, Zp,         // input
        xup1, xup2, xup3, 
        w1, w2, w3,  
        Phi1, Phi2, Phi3,   // output  
        Vp1, Vp2, Vp3,
        B1, B2, B3, 
        B4, B5, B6, 
        B7, B8, B9      
        );
}