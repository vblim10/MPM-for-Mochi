/*==========================================================
 * G2P_APIC_2d.c - grid to particle interpolation
 *
 * 1. Interpolates grid values (xup1,xup2) and (w1,w2)
 *          to particle values (Phi1,Phi2) and (V1,V2)
 *
 * 2. Updates APIC var Bp
 *
 * The calling syntax is:
 *
 *		[Phi1,Phi2, Vp1,Vp2, B1,B2,B3,B4] 
 *      = G2P_APIC_2d(Ox,Oy, hx,hy, nx,ny, Xp,Yp, xup1,xup2, w1,w2)
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

// G2P interpolation routine
void G2P(const double Ox, const double Oy,                  // scalars
         const double hx, const double hy, 
         size_t nx, size_t ny, size_t NP,
         double *Xp, double *Yp,                            // input
         double *xup1, double *xup2, 
         double *w1, double *w2,  
         double *Phi1, double *Phi2,                        // output
         double *Vp1, double *Vp2,     
         double *B1, double *B2, double *B3, double *B4
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
    double N[49];
    double bl1, bl2, br1, br2;
   
    // EVERYTHING CONDENSED IN THIS LOOP OVER PARTICLES
    for (int pidx=0; pidx<NP; pidx++) {
        // Reset temp vars
        fp_curr1 = 0.0; 
        fp_curr2 = 0.0;
        fp_curr3 = 0.0; 
        fp_curr4 = 0.0;

        // Position for particle p
        xp_curr = Xp[pidx];
        yp_curr = Yp[pidx];
      
        // Index of reference grid node for particle p 
        I = floor( (xp_curr-Ox) / hx );  
        J = floor( (yp_curr-Oy) / hy );  

        // 1st loop over stencil nodes
        for (int j=0; j<49; j++) {
            N[j] = 0.0;
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
                N[j] = N_weight(dx) * N_weight(dy);

                // PIC transfer
                fp_curr1 += xup1[gidx]*N[j];
                fp_curr2 += xup2[gidx]*N[j];
                fp_curr3 += w1[gidx]*N[j];
                fp_curr4 += w2[gidx]*N[j];
            }
        }
        Phi1[pidx] = fp_curr1;
        Phi2[pidx] = fp_curr2;
        Vp1[pidx] = fp_curr3;
        Vp2[pidx] = fp_curr4;
      
        // 2nd loop over stencil nodes (need a 2nd loop, bc we need Phi)
        for (int j=0; j<49; j++) {
            if (N[j] != 0) {
                k = floor(j/7);      // k: 0 to 6
                l = j - (7*k);       // l: 0 to 6
                k = I + k - 3;       // k: Ip + (-3 to 3)
                l = J + l - 3;       // l: Jp + (-3 to 3)
            
                // Compute gidx of (k,l) 
                gidx = (k)*(ny+1) + l; 
                xg_curr = Ox + hx*k;
                yg_curr = Oy + hy*l;

                bl1 = xg_curr - xp_curr + xup1[gidx] - Phi1[pidx];
                bl2 = yg_curr - yp_curr + xup2[gidx] - Phi2[pidx];
                br1 = xg_curr - xp_curr - xup1[gidx] + Phi1[pidx];
                br2 = yg_curr - yp_curr - xup2[gidx] + Phi2[pidx];

                // Compute APIC var Bp
                B1[pidx] += 0.5 * (w1[gidx] * bl1 +  br1 * w1[gidx]) *N[j];
                B2[pidx] += 0.5 * (w2[gidx] * bl1 +  br2 * w1[gidx]) *N[j];
                B3[pidx] += 0.5 * (w1[gidx] * bl2 +  br1 * w2[gidx]) *N[j];
                B4[pidx] += 0.5 * (w2[gidx] * bl2 +  br2 * w2[gidx]) *N[j];
            }
        }
    }
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    // Declare variables
    double *Xp, *Yp;                    // Particle positions (NP x 1)
    double *xup1, *xup2, *w1, *w2;      // grid values        (1 x NG)
    double *Phi1, *Phi2;                // particle values    (NP x 1)
    double *Vp1, *Vp2;                  //                    (1 x NP)
    double *B1, *B2, *B3, *B4;          // APIC var Bp        (1 x NP)
    //const double Ox; const double Oy;   // "origin" of Eul. grid
    //const double hx; const double hy;   // grid spacing
    size_t nx; size_t ny;               // Eul. grid partitions
    size_t NP;                          // # particles

    /* check for proper number of arguments */
    if (nrhs!=12) {
        mexErrMsgIdAndTxt("G2P_APIC_2d:nrhs","12 inputs required.");
    }
    if (nlhs!=8) {
        mexErrMsgIdAndTxt("G2P_APIC_2d:nlhs","8 outputs required.");
    }

    /* make sure the 6th input argument is type double */
    if (!mxIsDouble(prhs[6]) ||  mxIsComplex(prhs[6])) {
        mexErrMsgIdAndTxt("G2P_APIC_2d:notDouble","Xp must be type double.");
    }

    /* check that number of cols in input arguments is 1 */
    if (mxGetN(prhs[6])!=1 || mxGetN(prhs[7])!=1 ) {
        mexErrMsgIdAndTxt("G2P_APIC_2d:notRowVector","Xp, Yp must be column vectors.");
    }
   
    /* check that number of rows in input arguments is 1 */
    if (mxGetM(prhs[8])!=1 || mxGetM(prhs[9])!=1 ||
        mxGetM(prhs[10])!=1 || mxGetM(prhs[11])!=1  ) {
        mexErrMsgIdAndTxt("G2P_APIC_2d:notRowVector","xup1, xup2, w1, w2 must be row vectors.");
    }
   
   
   /* make sure the following input arguments are scalar */
    if (mxGetNumberOfElements(prhs[0])!=1 || mxGetNumberOfElements(prhs[1])!=1 ||
        mxGetNumberOfElements(prhs[2])!=1 || mxGetNumberOfElements(prhs[3])!=1 ||
        mxGetNumberOfElements(prhs[4])!=1 || mxGetNumberOfElements(prhs[5])!=1) {
        mexErrMsgIdAndTxt("G2P_APIC_2d:notScalar","Inputs Ox, Oy, hx, hy, nx, ny must be scalars.");
    }
   
    /* get the values of the scalar inputs  */
    const double Ox = mxGetScalar(prhs[0]);
    const double Oy = mxGetScalar(prhs[1]);
    const double hx = mxGetScalar(prhs[2]);
    const double hy = mxGetScalar(prhs[3]);
    nx = mxGetScalar(prhs[4]);
    ny = mxGetScalar(prhs[5]);
   
    /* create a pointer to the real data in Xp, Yp, xup1, xup2 */
    #if MX_HAS_INTERLEAVED_COMPLEX   // get a C error if we don't do this if/else statement
    Xp = mxGetDoubles(prhs[6]);
    Yp = mxGetDoubles(prhs[7]);
    xup1 = mxGetDoubles(prhs[8]);
    xup2 = mxGetDoubles(prhs[9]);
    w1 = mxGetDoubles(prhs[10]);
    w2 = mxGetDoubles(prhs[11]);
    #else
    Xp = mxGetPr(prhs[6]);           // matlab recommends not to use mxGetPr
    Yp = mxGetPr(prhs[7]);
    xup1 = mxGetPr(prhs[8]);
    xup2 = mxGetPr(prhs[9]);
    w1 = mxGetPr(prhs[10]);
    w2 = mxGetPr(prhs[11]);
    #endif
   
    /* get dimensions of the input matrix */
    NP = mxGetM(prhs[6]);

    /* create the output matrices */
    plhs[0] = mxCreateDoubleMatrix((mwSize)NP,1,mxREAL);
    plhs[1] = mxCreateDoubleMatrix((mwSize)NP,1,mxREAL);
    plhs[2] = mxCreateDoubleMatrix(1,(mwSize)NP,mxREAL);
    plhs[3] = mxCreateDoubleMatrix(1,(mwSize)NP,mxREAL);
    plhs[4] = mxCreateDoubleMatrix(1,(mwSize)NP,mxREAL);
    plhs[5] = mxCreateDoubleMatrix(1,(mwSize)NP,mxREAL);
    plhs[6] = mxCreateDoubleMatrix(1,(mwSize)NP,mxREAL);
    plhs[7] = mxCreateDoubleMatrix(1,(mwSize)NP,mxREAL);

    /* get a pointer to the real data in the output matrix */
    #if MX_HAS_INTERLEAVED_COMPLEX
    Phi1 = mxGetDoubles(plhs[0]);
    Phi2 = mxGetDoubles(plhs[1]);
    Vp1 = mxGetDoubles(plhs[2]);
    Vp2 = mxGetDoubles(plhs[3]);
    B1 = mxGetDoubles(plhs[4]);
    B2 = mxGetDoubles(plhs[5]);
    B3 = mxGetDoubles(plhs[6]);
    B4 = mxGetDoubles(plhs[7]);
    #else
    Phi1 = mxGetPr(plhs[0]); 
    Phi2 = mxGetPr(plhs[1]);
    Vp1 = mxGetPr(plhs[2]); 
    Vp2 = mxGetPr(plhs[3]);
    B1 = mxGetPr(plhs[4]);
    B2 = mxGetPr(plhs[5]);
    B3 = mxGetPr(plhs[6]);
    B4 = mxGetPr(plhs[7]);
    #endif
   
    /* call the computational routine */
    G2P(Ox, Oy, hx, hy, nx, ny, NP,     // scalars
        Xp, Yp,                         // input
        xup1, xup2, 
        w1, w2,     
        Phi1, Phi2,                     // output
        Vp1, Vp2,
        B1, B2, B3, B4
        );
}