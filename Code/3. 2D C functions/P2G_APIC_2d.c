/*==========================================================
 * P2G_APIC_2d.c - grid to particle interpolation
 *
 * Interpolates particle values mom_p = [mom_p1; mom_p2] & mass_p
 *               to grid values mom_g = [mom_g1; mom_g2] & mass_g
 *
 * The calling syntax is:
 *
 *		[mom_g1,mom_g2, mass_g] 
 *      = P2G_APIC_2d(Ox,Oy, hx,hy, nx,ny, Xp,Yp, B1,B2,B3,B4, mom_p1,mom_p2, mass_p)
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

// P2G interpolation routine
void P2G(const double Ox, const double Oy,                 // scalars
         const double hx, const double hy, 
         size_t nx, size_t ny, size_t NP,
         double *Xp, double *Yp,                           // input
         double *B1, double *B2, double *B3, double *B4,
         double *mom_p1, double *mom_p2, double *mass_p,
         double *mom_g1, double *mom_g2, double *mass_g    // output      
         )
{
    // Constants
    const int d = 2;

    // Temp vars
    double xp_curr, yp_curr;
    double xg_curr, yg_curr, dx, dy;
    int I, J;
    int l, k, gidx;
    double N[49];
    double D1, D2, D3, D4;        // Dp = [D1 D3; D2 D4]
    double C1, C2, C3, C4, c;     // Cp = [C1 C3; C2 C4]
    double A1, A2;                // [A1; A2] = Cp * ["dx"; "dy"]
   
    // EVERYTHING CONDENSED IN THIS LOOP OVER PARTICLES
    for (int pidx=0; pidx<NP; pidx++) {
        // Reset temp vars
        D1=0.0; D2=0.0; D3=0.0; D4=0.0; 

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
                mass_g[gidx] += mass_p[pidx]*N[j];
            
                // Compute APIC var Dp 
                dx = xg_curr - xp_curr; // *** not same dx,dy used to compute N ***
                dy = yg_curr - yp_curr;  
                D1 += (dx * dx * N[j]);
                D2 += (dx * dy * N[j]);
                D3 += (dx * dy * N[j]);
                D4 += (dy * dy * N[j]);
            }
        }
      
        // APIC var Cp 
        c = 1.0/(D1*D4 - D2*D3);
        C1 = c * ( B1[pidx]*D4 - B3[pidx]*D2);
        C2 = c * ( B2[pidx]*D4 - B4[pidx]*D2);
        C3 = c * (-B1[pidx]*D3 + B3[pidx]*D1);
        C4 = c * (-B2[pidx]*D3 + B4[pidx]*D1);
      
        // 2nd loop over stencil nodes (since we need Cp)
        for (int j=0; j<49; j++) {
            if (N[j] != 0) {
                // can just make array of gidx,dx,dy values like N[j]
                // but simple arith. is fast enough (& less storage)?
                k = floor(j/7);      // k: 0 to 6
                l = j - (7*k);       // l: 0 to 6
                k = I + k - 3;       // k: Ip + (-3 to 3)
                l = J + l - 3;       // l: Jp + (-3 to 3)

                // Compute gidx of (k,l) 
                gidx = (k)*(ny+1) + l; 
                xg_curr = Ox + hx*k;
                yg_curr = Oy + hy*l;
                dx = xg_curr - xp_curr; // *** not same dx,dy used to compute N ***
                dy = yg_curr - yp_curr; 

                // APIC var Cp*X "angular velocity"
                A1 = C1*dx + C3*dy;
                A2 = C2*dx + C4*dy;

                // APIC transfer
                mom_g1[gidx] += ( mom_p1[pidx] + (mass_p[pidx] * A1) ) * N[j];
                mom_g2[gidx] += ( mom_p2[pidx] + (mass_p[pidx] * A2) ) * N[j];
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
    double *B1, *B2, *B3, *B4;          // APIC vars          (1 x NP), where Bp = [B1;B2;B3;B4] 
    double *mom_p1, *mom_p2, *mass_p;   // particle values    (1 x NP), (1 x NP), (NP x 1)
    double *mom_g1, *mom_g2, *mass_g;   // grid values        (1 x NG)
    //const double Ox; const double Oy;   // "origin" of Eul. grid
    //const double hx; const double hy;   // grid spacing
    size_t nx; size_t ny;               // Eul. grid partitions
    size_t NP;                          // # particles

    /* check for proper number of arguments */
    if (nrhs!=15) {
        mexErrMsgIdAndTxt("P2G_APIC_2d:nrhs","15 inputs required.");
    }
    if (nlhs!=3) {
        mexErrMsgIdAndTxt("P2G_APIC_2d:nlhs","3 outputs required.");
    }

    /* make sure the 6th input argument is type double */
    if (!mxIsDouble(prhs[6]) ||  mxIsComplex(prhs[6])) {
        mexErrMsgIdAndTxt("P2G_APIC_2d:notDouble","Xp must be type double.");
    }

    /* check that number of cols in input arguments is 1 */
    if (mxGetN(prhs[6])!=1 || mxGetN(prhs[7])!=1 || mxGetN(prhs[14])!=1) {
        mexErrMsgIdAndTxt("P2G_APIC_2d:notRowVector","Xp, Yp, mass_p must be a column vectors.");
    }
   
    /* check that number of rows in input arguments is 1 */
    if (mxGetM(prhs[8])!=1 || mxGetM(prhs[9])!=1 || mxGetM(prhs[10])!=1 || mxGetM(prhs[11])!=1 ||
        mxGetM(prhs[12])!=1 || mxGetM(prhs[13])!=1) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:notRowVector","B1, B2, B3, B4, mom_p1, mom_p2 must be a row vectors.");
    }
   
    /* make sure the following input arguments are scalar */
    if (mxGetNumberOfElements(prhs[0])!=1 || mxGetNumberOfElements(prhs[1])!=1 ||
        mxGetNumberOfElements(prhs[2])!=1 || mxGetNumberOfElements(prhs[3])!=1 ||
        mxGetNumberOfElements(prhs[4])!=1 || mxGetNumberOfElements(prhs[5])!=1) {
        mexErrMsgIdAndTxt("P2G_APIC_2d:notScalar","Inputs Ox, Oy, hx, hy, nx, ny must be scalars.");
    }
   
    /* get the values of the scalar inputs  */
    const double Ox = mxGetScalar(prhs[0]);
    const double Oy = mxGetScalar(prhs[1]);
    const double hx = mxGetScalar(prhs[2]);
    const double hy = mxGetScalar(prhs[3]);
    nx = mxGetScalar(prhs[4]);
    ny = mxGetScalar(prhs[5]);
   
    /* create a pointer to the real data in Xp, Yp, mom_p1 */
    #if MX_HAS_INTERLEAVED_COMPLEX      // get a C error if we don't do this if/else statement
    Xp = mxGetDoubles(prhs[6]);
    Yp = mxGetDoubles(prhs[7]);
    B1 = mxGetDoubles(prhs[8]);
    B2 = mxGetDoubles(prhs[9]);
    B3 = mxGetDoubles(prhs[10]);
    B4 = mxGetDoubles(prhs[11]);
    mom_p1 = mxGetDoubles(prhs[12]);
    mom_p2 = mxGetDoubles(prhs[13]);
    mass_p = mxGetDoubles(prhs[14]);
    #else
    Xp = mxGetPr(prhs[6]);        // matlab recommends not to use mxGetPr
    Yp = mxGetPr(prhs[7]);
    B1 = mxGetPr(prhs[8]);
    B2 = mxGetPr(prhs[9]);
    B3 = mxGetPr(prhs[10]);
    B4 = mxGetPr(prhs[11]);
    mom_p1 = mxGetPr(prhs[12]);
    mom_p2 = mxGetPr(prhs[13]);
    mass_p = mxGetPr(prhs[14]);
    #endif
   
    /* get dimensions of the input matrix */
    NP = mxGetM(prhs[6]);
    size_t NG = (nx+1)*(ny+1);

    /* create the output matrix */
    plhs[0] = mxCreateDoubleMatrix(1,(mwSize)NG,mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1,(mwSize)NG,mxREAL);
    plhs[2] = mxCreateDoubleMatrix((mwSize)NG,1,mxREAL);

    /* get a pointer to the real data in the output matrix */
    #if MX_HAS_INTERLEAVED_COMPLEX
    mom_g1 = mxGetDoubles(plhs[0]);
    mom_g2 = mxGetDoubles(plhs[1]);
    mass_g = mxGetDoubles(plhs[2]);
    #else
    mom_g1 = mxGetPr(plhs[0]);
    mom_g2 = mxGetPr(plhs[1]);
    mass_g = mxGetPr(plhs[2]);
    #endif
   
    /* call the computational routine */
    P2G(Ox, Oy, hx, hy, nx, ny, NP, // scalars
        Xp, Yp,                     // input
        B1, B2, B3, B4,
        mom_p1, mom_p2, mass_p,
        mom_g1, mom_g2, mass_g      // output 
        );
}