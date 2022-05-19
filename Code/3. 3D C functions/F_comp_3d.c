/*==========================================================
 * F_comp_3d.c - grid to particle interpolation
 *
 * Computes deformation gradient F 
 *
 * The calling syntax is:
 *
 *		[Fup1, Fup2, Fup3, Fup4, Fup5, Fup6, Fup7, Fup8, Fup9] 
 *      = F_comp_2d(Ox,Oy,Oz, hx,hy,hz, nx,ny,nz, Xp,Yp,Zp, dt, vx,vy,vz, F1,F2,F3,F4,F5,F6,F7,F8,F9)
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
void F_comp(const double Ox, const double Oy, const double Oz, // constants
            const double hx, const double hy, const double hz,
            size_t nx, size_t ny, size_t nz, size_t NP,
            double *Xp, double *Yp, double *Zp,                // input
            double dt, double *vx, double *vy, double *vz,
            double *F1, double *F2, double *F3, 
            double *F4, double *F5, double *F6,
            double *F7, double *F8, double *F9,
            double *Fup1, double *Fup2, double *Fup3,          // output 
            double *Fup4, double *Fup5, double *Fup6,
            double *Fup7, double *Fup8, double *Fup9
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
        // Ft = Id + dt*Ft
        fp1_curr = dt*fp1_curr + 1;
        fp2_curr = dt*fp2_curr;
        fp3_curr = dt*fp3_curr;
        
        fp4_curr = dt*fp4_curr;
        fp5_curr = dt*fp5_curr + 1;
        fp6_curr = dt*fp6_curr;
        
        fp7_curr = dt*fp7_curr;
        fp8_curr = dt*fp8_curr;
        fp9_curr = dt*fp9_curr + 1;
        
        // Fup = Ft * Fkp
        Fup1[pidx] = fp1_curr*F1[pidx] + fp4_curr*F2[pidx] + fp7_curr*F3[pidx];
        Fup2[pidx] = fp2_curr*F1[pidx] + fp5_curr*F2[pidx] + fp8_curr*F3[pidx];
        Fup3[pidx] = fp3_curr*F1[pidx] + fp6_curr*F2[pidx] + fp9_curr*F3[pidx];

        Fup4[pidx] = fp1_curr*F4[pidx] + fp4_curr*F5[pidx] + fp7_curr*F6[pidx];
        Fup5[pidx] = fp2_curr*F4[pidx] + fp5_curr*F5[pidx] + fp8_curr*F6[pidx];
        Fup6[pidx] = fp3_curr*F4[pidx] + fp6_curr*F5[pidx] + fp9_curr*F6[pidx];

        Fup7[pidx] = fp1_curr*F7[pidx] + fp4_curr*F8[pidx] + fp7_curr*F9[pidx];
        Fup8[pidx] = fp2_curr*F7[pidx] + fp5_curr*F8[pidx] + fp8_curr*F9[pidx];
        Fup9[pidx] = fp3_curr*F7[pidx] + fp6_curr*F8[pidx] + fp9_curr*F9[pidx];
    }
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    // Declare variables
    double *Xp; double *Yp; double *Zp;     // Particle positions (NP x 1)
    double dt;                              // time step
    double *vx,   *vy,   *vz;               // grid vel v=(u+w)   (NG x 1)
    double *F1,   *F4,   *F7;                   
    double *F2,   *F5,   *F8;               // def. grad. Fpk     (NP x 1)
    double *F3,   *F6,   *F9;
    double *Fup1, *Fup4, *Fup7;                   
    double *Fup2, *Fup5, *Fup8;             // F update           (NP x 1)
    double *Fup3, *Fup6, *Fup9;
    //const double Ox; const double Oy; const double Oz;  // "origin" of Eul. grid
    //const double hx; const double hy; const double hz;  // grid spacing
    size_t nx; size_t ny; size_t nz;        // Eul. grid partitions
    size_t NP;                              // # particles

    /* check for proper number of arguments */
    if (nrhs!=25) {
        mexErrMsgIdAndTxt("F_comp_3d:nrhs","25 inputs required.");
    }
    if (nlhs!=9) {
        mexErrMsgIdAndTxt("F_comp_3d:nlhs","9 outputs required.");
    }

    /* make sure the 10th input argument is type double */
    if ( !mxIsDouble(prhs[9]) ||  mxIsComplex(prhs[9])) {
        mexErrMsgIdAndTxt("F_comp_3d:notDouble","Xp must be type double.");
    }

    /* check that number of cols in input arguments is 1 */
    if (mxGetN(prhs[9])!=1 || mxGetN(prhs[10])!=1 || mxGetN(prhs[11])!=1) {
        mexErrMsgIdAndTxt("F_comp_3d:notColVector","Xp, Yp, Zp must be column vectors.");
    }
   
    /* check that number of rows in input arguments is 1 */
    if (mxGetM(prhs[16])!=1 || mxGetM(prhs[17])!=1 || mxGetM(prhs[18])!=1 || 
        mxGetM(prhs[19])!=1 || mxGetM(prhs[20])!=1 || mxGetM(prhs[21])!=1 || 
        mxGetM(prhs[22])!=1 || mxGetM(prhs[23])!=1 || mxGetM(prhs[24])!=1 ) {
        mexErrMsgIdAndTxt("F_comp_3d:notRowVector","F1, F2, ... , F8, F9 must be row vectors.");
    }
   
    /* make sure the following input arguments are scalar */
    if (mxGetNumberOfElements(prhs[0])!=1  || mxGetNumberOfElements(prhs[1])!=1  ||
        mxGetNumberOfElements(prhs[2])!=1  || mxGetNumberOfElements(prhs[3])!=1  ||
        mxGetNumberOfElements(prhs[4])!=1  || mxGetNumberOfElements(prhs[5])!=1  ||
        mxGetNumberOfElements(prhs[6])!=1  || mxGetNumberOfElements(prhs[7])!=1  ||
        mxGetNumberOfElements(prhs[8])!=1  || mxGetNumberOfElements(prhs[12])!=1 ) {
        mexErrMsgIdAndTxt("F_comp_3d:notScalar","Inputs Ox, Oy, Oz, hx, hy, hz, nx, ny, nz, dt must be scalars.");
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
    dt = mxGetScalar(prhs[12]);
   
   
    /* create a pointer to the real data in Xp, Yp, Zp, F1, F2, ..., F8, F9 */
    #if MX_HAS_INTERLEAVED_COMPLEX      // get a C error if we don't do this if/else statement
    Xp = mxGetDoubles(prhs[9]);
    Yp = mxGetDoubles(prhs[10]);
    Zp = mxGetDoubles(prhs[11]);
    vx = mxGetDoubles(prhs[13]);
    vy = mxGetDoubles(prhs[14]);
    vz = mxGetDoubles(prhs[15]);
    F1 = mxGetDoubles(prhs[16]);
    F2 = mxGetDoubles(prhs[17]);
    F3 = mxGetDoubles(prhs[18]);
    F4 = mxGetDoubles(prhs[19]);
    F5 = mxGetDoubles(prhs[20]);
    F6 = mxGetDoubles(prhs[21]);
    F7 = mxGetDoubles(prhs[22]);
    F8 = mxGetDoubles(prhs[23]);
    F9 = mxGetDoubles(prhs[24]);
    #else
    Xp = mxGetPr(prhs[9]);
    Yp = mxGetPr(prhs[10]);
    Zp = mxGetPr(prhs[11]);
    vx = mxGetPr(prhs[13]);
    vy = mxGetPr(prhs[14]);
    vz = mxGetPr(prhs[15]);
    F1 = mxGetPr(prhs[16]);
    F2 = mxGetPr(prhs[17]);
    F3 = mxGetPr(prhs[18]);
    F4 = mxGetPr(prhs[19]);
    F5 = mxGetPr(prhs[20]);
    F6 = mxGetPr(prhs[21]);
    F7 = mxGetPr(prhs[22]);
    F8 = mxGetPr(prhs[23]);
    F9 = mxGetPr(prhs[24]);
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
    Fup1 = mxGetDoubles(plhs[0]);
    Fup2 = mxGetDoubles(plhs[1]);
    Fup3 = mxGetDoubles(plhs[2]);
    Fup4 = mxGetDoubles(plhs[3]);
    Fup5 = mxGetDoubles(plhs[4]);
    Fup6 = mxGetDoubles(plhs[5]);
    Fup7 = mxGetDoubles(plhs[6]);
    Fup8 = mxGetDoubles(plhs[7]);
    Fup9 = mxGetDoubles(plhs[8]);
    #else
    Fup1 = mxGetPr(plhs[0]);
    Fup2 = mxGetPr(plhs[1]);
    Fup3 = mxGetPr(plhs[2]);
    Fup4 = mxGetPr(plhs[3]);
    Fup5 = mxGetPr(plhs[4]);
    Fup6 = mxGetPr(plhs[5]);
    Fup7 = mxGetPr(plhs[6]);
    Fup8 = mxGetPr(plhs[7]);
    Fup9 = mxGetPr(plhs[8]);
    #endif
   
    /* call the computational routine */
    F_comp(Ox, Oy, Oz,          // constants
           hx, hy, hz,
           nx, ny, nz, NP,
           Xp, Yp, Zp,          // input
           dt, vx, vy, vz,
           F1, F2, F3, 
           F4, F5, F6,
           F7, F8, F9,
           Fup1, Fup2, Fup3,    // output 
           Fup4, Fup5, Fup6,
           Fup7, Fup8, Fup9
           );
}