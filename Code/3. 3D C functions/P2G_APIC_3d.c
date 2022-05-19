/*==========================================================
 * P2G_APIC_3d.c - grid to particle interpolation
 *
 * Interpolates particle values mom_p = [mom_p1; mom_p2; mom_p3] & mass_p
 *               to grid values mom_g = [mom_g1; mom_g2; mom_g3] & mass_g
 *
 * The calling syntax is:
 *
 *		[mom_g1,mom_g2,mom_g3, mass_g] 
 *      = P2G_APIC_3d(Ox,Oy,Oz, hx,hy,hz, nx,ny,nz, Xp,Yp,Zp, 
 *                    B1,B2,B3,B4,B5,B6,B7,B8,B9, mom_p1,mom_p2,mom_p3, mass_p)
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

// P2G interpolation routine
void P2G(const double Ox, const double Oy, const double Oz, // constants
         const double hx, const double hy, const double hz,
         size_t nx, size_t ny, size_t nz, size_t NP,
         double *Xp, double *Yp, double *Zp,                // input
         double *B1, double *B2, double *B3, 
         double *B4, double *B5, double *B6, 
         double *B7, double *B8, double *B9,
         double *mom_p1, double *mom_p2, double *mom_p3, 
         double *mass_p,
         double *mom_g1, double *mom_g2, double *mom_g3,    // output  
         double *mass_g       
        )
{
    // Constants
    const int d = 3;
   
    // Temp vars
    double xp_curr, yp_curr, zp_curr;
    double xg_curr, yg_curr, zg_curr, dx, dy, dz;
    int I, J, K;
    int i, j, k, gidx;
    double D1, D2, D3;  //      [ D1 D4 D7 ]
    double D4, D5, D6;  // Dp = | D2 D5 D8 |
    double D7, D8, D9;  //      [ D3 D6 D9 ]
    double detD; 
    double E1, E2, E3;  
    double E4, E5, E6;  // E = inv(Dp)
    double E7, E8, E9;  
    double C1, C2, C3;  //      [ C1 C4 C7 ]
    double C4, C5, C6;  // Cp = | C2 C5 C8 |
    double C7, C8, C9;  //      [ C3 C6 C9 ]
    double a1, a2, a3;  // a = [a1; a2; a3] = Cp * ["dx"; "dy"; "dz"]
   
    // EVERYTHING CONDENSED IN THIS LOOP OVER PARTICLES
    for (int pidx=0; pidx<NP; pidx++) {
        // Reset temp vars
        double N[7][49] = {0.0}; // array of all zeros
        D1=0.0; D2=0.0; D3=0.0; 
        D4=0.0; D5=0.0; D6=0.0;
        D7=0.0; D8=0.0; D9=0.0;
        
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
                        mass_g[gidx] += mass_p[pidx]*N[level][S];
                        
                        // Compute APIC var Dp 
                        dx = xg_curr - xp_curr; // *** not same dx,dy,dz used to compute N ***
                        dy = yg_curr - yp_curr;  
                        dz = zg_curr - zp_curr;
                        
                        D1 += (dx * dx * N[level][S]);
                        D2 += (dy * dx * N[level][S]);
                        D3 += (dz * dx * N[level][S]);
                        
                        D4 += (dx * dy * N[level][S]);
                        D5 += (dy * dy * N[level][S]);
                        D6 += (dz * dy * N[level][S]);
                        
                        D7 += (dx * dz * N[level][S]);
                        D8 += (dy * dz * N[level][S]);
                        D9 += (dz * dz * N[level][S]);
                    }
                }
            }
        }
        // Compute APIC var Cp = Bp*E
        detD = D1*(D5*D9-D8*D6) - D4*(D2*D9-D8*D3) + D7*(D2*D6-D5*D3);
        E1 = (D5*D9 - D8*D6) / detD;
        E2 = (D8*D3 - D2*D9) / detD;
        E3 = (D2*D6 - D5*D3) / detD;
        
        E4 = (D7*D6 - D4*D9) / detD;
        E5 = (D1*D9 - D7*D3) / detD;
        E6 = (D4*D3 - D1*D6) / detD;
        
        E7 = (D4*D8 - D7*D5) / detD;
        E8 = (D7*D2 - D1*D8) / detD;
        E9 = (D1*D5 - D4*D2) / detD;
        
        C1 = B1[pidx]*E1 + B4[pidx]*E4 + B7[pidx]*E7;
        C2 = B2[pidx]*E1 + B5[pidx]*E4 + B8[pidx]*E7;
        C3 = B3[pidx]*E1 + B6[pidx]*E4 + B9[pidx]*E7;

        C4 = B1[pidx]*E2 + B4[pidx]*E5 + B7[pidx]*E8;
        C5 = B2[pidx]*E2 + B5[pidx]*E5 + B8[pidx]*E8;
        C6 = B3[pidx]*E2 + B6[pidx]*E5 + B9[pidx]*E8;

        C7 = B1[pidx]*E3 + B4[pidx]*E6 + B7[pidx]*E9;
        C8 = B2[pidx]*E3 + B5[pidx]*E6 + B8[pidx]*E9;
        C9 = B3[pidx]*E3 + B6[pidx]*E6 + B9[pidx]*E9;
        
        // 2nd (*) loop (need a 2nd double-loop b/c we need Cp)
        // Loop over each level k
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
                    
                    // APIC var a = Cp * ["dx"; "dy"; "dz"] "angular velocity"
                    dx = xg_curr - xp_curr; // *** not same dx,dy,dz used to compute N ***
                    dy = yg_curr - yp_curr; 
                    dz = zg_curr - zp_curr;
                    a1 = C1*dx + C4*dy + C7*dz;
                    a2 = C2*dx + C5*dy + C8*dz; 
                    a3 = C3*dx + C6*dy + C9*dz;
                    
                    // APIC transfer
                    mom_g1[gidx] += ( mom_p1[pidx] + (mass_p[pidx] * a1) ) * N[level][S];
                    mom_g2[gidx] += ( mom_p2[pidx] + (mass_p[pidx] * a2) ) * N[level][S];
                    mom_g3[gidx] += ( mom_p3[pidx] + (mass_p[pidx] * a3) ) * N[level][S];
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
    double *B1, *B2, *B3;               // APIC vars          (1 x NP) 
    double *B4, *B5, *B6;               // Bp = [B1 B4 B7; B2 B5 B8; B3 B6 B9]
    double *B7, *B8, *B9;               
    double *mom_p1, *mom_p2, *mom_p3;   // particle momentum  (1 x NP)
    double *mass_p;                     // particle mass      (NP x 1)
    double *mom_g1, *mom_g2, *mom_g3;   // grid momentum      (1 x NG)
    double *mass_g;                     // grid mass          (NG x 1)
    //const double Ox; const double Oy; const double Oz;  // "origin" of Eul. grid
    //const double hx; const double hy; const double hz;  // grid spacing
    size_t nx; size_t ny; size_t nz;    // Eul. grid partitions
    size_t NP;                          // # particles

    /* check for proper number of arguments */
    if (nrhs!=25) {
        mexErrMsgIdAndTxt("P2G_APIC_3d:nrhs","25 inputs required.");
    }
    if (nlhs!=4) {
        mexErrMsgIdAndTxt("P2G_APIC_3d:nlhs","4 outputs required.");
    }

    /* make sure the 10th input argument is type double */
    if (!mxIsDouble(prhs[9]) ||  mxIsComplex(prhs[9])) {
        mexErrMsgIdAndTxt("P2G_APIC_3d:notDouble","Xp must be type double.");
    }

    /* check that number of cols in input arguments is 1 */
    if (mxGetN(prhs[9])!=1 || mxGetN(prhs[10])!=1 || mxGetN(prhs[11])!=1 || mxGetN(prhs[24])!=1) {
        mexErrMsgIdAndTxt("P2G_APIC_3d:notColVector","Xp, Yp, Zp, mass_p must be column vectors.");
    }
   
    /* check that number of rows in input arguments is 1 */
    if (mxGetM(prhs[12])!=1 || mxGetM(prhs[13])!=1 || mxGetM(prhs[14])!=1 || 
        mxGetM(prhs[15])!=1 || mxGetM(prhs[16])!=1 || mxGetM(prhs[17])!=1 || 
        mxGetM(prhs[18])!=1 || mxGetM(prhs[19])!=1 || mxGetM(prhs[20])!=1 || 
        mxGetM(prhs[21])!=1 || mxGetM(prhs[22])!=1 || mxGetM(prhs[23])!=1 ) {
        mexErrMsgIdAndTxt("P2G_APIC_3d:notRowVector","B1, B2,..., B8, B9, mom_p1, mom_p2, mom_p3 must be row vectors.");
    }
   
    /* make sure the following input arguments are scalar */
    if (mxGetNumberOfElements(prhs[0])!=1  || mxGetNumberOfElements(prhs[1])!=1  ||
        mxGetNumberOfElements(prhs[2])!=1  || mxGetNumberOfElements(prhs[3])!=1  ||
        mxGetNumberOfElements(prhs[4])!=1  || mxGetNumberOfElements(prhs[5])!=1  ||
        mxGetNumberOfElements(prhs[6])!=1  || mxGetNumberOfElements(prhs[7])!=1  ||
        mxGetNumberOfElements(prhs[8])!=1  ) {
        mexErrMsgIdAndTxt("P2G_APIC_3d:notScalar","Inputs Ox, Oy, Oz, hx, hy, hz, nx, ny, nz must be scalars.");
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
   
    /* create a pointer to the real data in Xp, Yp, Zp, fp1, fp2, fp3 */
    #if MX_HAS_INTERLEAVED_COMPLEX      // get a C error if we don't do this if/else statement
    Xp = mxGetDoubles(prhs[9]);
    Yp = mxGetDoubles(prhs[10]);
    Zp = mxGetDoubles(prhs[11]);
    B1 = mxGetDoubles(prhs[12]);
    B2 = mxGetDoubles(prhs[13]);
    B3 = mxGetDoubles(prhs[14]);
    B4 = mxGetDoubles(prhs[15]);
    B5 = mxGetDoubles(prhs[16]);
    B6 = mxGetDoubles(prhs[17]);
    B7 = mxGetDoubles(prhs[18]);
    B8 = mxGetDoubles(prhs[19]);
    B9 = mxGetDoubles(prhs[20]);
    mom_p1 = mxGetDoubles(prhs[21]);
    mom_p2 = mxGetDoubles(prhs[22]);
    mom_p3 = mxGetDoubles(prhs[23]);
    mass_p = mxGetDoubles(prhs[24]);
    #else
    Xp = mxGetPr(prhs[9]);              // matlab recommends not to use mxGetPr
    Yp = mxGetPr(prhs[10]);
    Zp = mxGetPr(prhs[11]);
    B1 = mxGetPr(prhs[12]);
    B2 = mxGetPr(prhs[13]);
    B3 = mxGetPr(prhs[14]);
    B4 = mxGetPr(prhs[15]);
    B5 = mxGetPr(prhs[16]);
    B6 = mxGetPr(prhs[17]);
    B7 = mxGetPr(prhs[18]);
    B8 = mxGetPr(prhs[19]);
    B9 = mxGetPr(prhs[20]);
    mom_p1 = mxGetPr(prhs[21]);
    mom_p2 = mxGetPr(prhs[22]);
    mom_p3 = mxGetPr(prhs[23]);
    mass_p = mxGetPr(prhs[24]);
    #endif
   
    /* get dimensions of the input/output matrices */
    NP = mxGetM(prhs[9]);
    size_t NG = (nx+1)*(ny+1)*(nz+1);

    /* create the output matrices */
    plhs[0] = mxCreateDoubleMatrix(1,(mwSize)NG,mxREAL);
    plhs[1] = mxCreateDoubleMatrix(1,(mwSize)NG,mxREAL);
    plhs[2] = mxCreateDoubleMatrix(1,(mwSize)NG,mxREAL);
    plhs[3] = mxCreateDoubleMatrix((mwSize)NG,1,mxREAL);

    /* get a pointer to the real data in the output matrix */
    #if MX_HAS_INTERLEAVED_COMPLEX
    mom_g1 = mxGetDoubles(plhs[0]);
    mom_g2 = mxGetDoubles(plhs[1]);
    mom_g3 = mxGetDoubles(plhs[2]);
    mass_g = mxGetDoubles(plhs[3]);
    #else
    mom_g1 = mxGetPr(plhs[0]);
    mom_g2 = mxGetPr(plhs[1]);
    mom_g3 = mxGetPr(plhs[2]);
    mass_g = mxGetPr(plhs[3]);
    #endif
   
    /* call the computational routine */
    P2G(Ox, Oy, Oz,             // constants
        hx, hy, hz,
        nx, ny, nz, NP,
        Xp, Yp, Zp,             // input
        B1, B2, B3, 
        B4, B5, B6, 
        B7, B8, B9,
        mom_p1, mom_p2, mom_p3, 
        mass_p,
        mom_g1, mom_g2, mom_g3, // output  
        mass_g       
        );
}