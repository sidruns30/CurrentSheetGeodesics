#include "metric.hpp"

// 0 for Schwarzschild, 1 for Boyer-Lundquist,  2 for Kerr Schild and 3 for MKS
void GcovFunc(int metric_type, double X[NDIM], double gcov[NDIM][NDIM])
{
    MUNULOOP gcov[mu][nu] = 0.;
    double r = X[1];
    double theta = X[2];

    if (metric_type == 0)
    {
        GetSCCovariantMetric(gcov, X);
    }

    else if (metric_type == 1)
    {
        GetBLCovariantMetric(gcov, X);
    }

    else if (metric_type == 2)
    {
        GetKSCovariantMetric(gcov, X);
    }

    else if (metric_type == 3)
    {
		// Make sure that the coordinates X are Kerr Schild
        GetMKSCovariantMetric_v2(gcov, X);
    }

    return;
}


void GconFunc(int metric_type, double X[NDIM], double gcon[NDIM][NDIM])
{
    MUNULOOP gcon[mu][nu] = 0.;
    double r = X[1];
    double theta = X[2];

    if (metric_type == 0)
    {
        GetSCContravariantMetric(gcon, X);
    }

    else if (metric_type == 1)
    {
        GetBLContravariantMetric(gcon, X);
    }

    else if (metric_type == 2)
    {
        GetKSContravariantMetric(gcon, X);
    }

	else if (metric_type == 3)
	{
		// Make sure that the coordinates are Kerr Schild
		GetMKSContravariantMetric_v2(gcon, X);
	}
    return;
}


#define DEL (1.e-7)
void GetConnectionNumerical(int metric_type, double X[NDIM], double conn[NDIM][NDIM][NDIM])
{
    double tmp[NDIM][NDIM][NDIM];
    double Xh[NDIM], Xl[NDIM];
    double gcon[NDIM][NDIM];
    double gcov[NDIM][NDIM];
    double gh[NDIM][NDIM];
    double gl[NDIM][NDIM];

    GcovFunc(metric_type, X, gcov);
    GconFunc(metric_type, X, gcon);

    for (int k = 0; k < NDIM; k++) 
    {
        for (int l = 0; l < NDIM; l++)
        {
            Xh[l] = X[l];
        }

        for (int l = 0; l < NDIM; l++)
        {
            Xl[l] = X[l];
        }

        Xh[k] += DEL;
        Xl[k] -= DEL;
        GcovFunc(metric_type, Xh, gh);
        GcovFunc(metric_type, Xl, gl);

        for (int i = 0; i < NDIM; i++) 
        {
            for (int j = 0; j < NDIM; j++)
            {
                conn[i][j][k] = (gh[i][j] - gl[i][j]) / (Xh[k] - Xl[k]);
            }
        }
    }

    // Rearrange to find \Gamma_{ijk}
    for (int i = 0; i < NDIM; i++)
    {
        for (int j = 0; j < NDIM; j++)
        {
            for (int k = 0; k < NDIM; k++)
            {
                tmp[i][j][k] = 0.5 * (conn[j][i][k] + conn[k][i][j] - conn[k][j][i]);
            }
        }
    }

    // G_{ijk} -> G^i_{jk}
    for (int i = 0; i < NDIM; i++) 
    {
        for (int j = 0; j < NDIM; j++) 
        {
            for (int k = 0; k < NDIM; k++) 
            {
                conn[i][j][k] = 0.;
                for (int l = 0; l < NDIM; l++)
                {
                    conn[i][j][k] += gcon[i][l] * tmp[l][j][k];
                }
            }
        }
    }

    return;
}
#undef DEL

// If the metric is known
void GetConnectionAnalytic(int metric_type, double X[NDIM], double gamma[NDIM][NDIM][NDIM])
{
    int i, j, k;
    for (i=0;i<NDIM;i++)
    {
        for (j=0;j<NDIM;j++)
        {
            for (k=0; k<NDIM; k++)
            {
                gamma[i][j][k] = 0.;
            }
        }
    }
    // Schwarzschild
    if (metric_type == 0)
    {
        GetSCConnection(X, gamma);
    }

    // Boyer-Lindquist
    else if (metric_type == 1)
    {
        GetBLConnection(X, gamma);
    }

    // Kerr Schild
    else if (metric_type == 2)
    {
        GetKSConnection(X, gamma);
    }

    // Modified Kerr Schild
    else if (metric_type == 3)
    {
        GetMKSConnection(X, gamma);
    }

    else
    {
        throw std::invalid_argument("metric type must be between [0-3]");
    }
}

// Go from metric type 1 to 2. Metric types [1, 2, 3]: BL, KS, MKS
void TransformCoordinates(int metric_type_1, int metric_type_2, double X1[NDIM], double X2[NDIM])
{
    // BL to KS / MKS
    if (metric_type_1 == 1)
    {
        // KS
        if (metric_type_2 == 2)
        {
            X_BLToKS(X1, X2);
            return;
        }

        // MKS
        else if (metric_type_2 == 3)
        {
            double XKS[NDIM];
            X_BLToKS(X1, XKS);
            X_KSToMKS_v2(XKS, X2);
            return;
        }

        else
        {
            throw std::invalid_argument("Only BL to KS and MKS supported");
            return;
        }
    }

    // KS to BL / MKS
    else if (metric_type_1 == 2)
    {
        // BL
        if (metric_type_2 == 1)
        {
            X_KSToBL(X1, X2);
            return;
        }

        // MKS
        else if (metric_type_2 == 3)
        {
            X_KSToMKS_v2(X1, X2);
            return;
        }

        else
        {
            throw std::invalid_argument("Only KS to BL and MKS supported");
            return;
        }
    }

    // MKS to BL / KS
    else if (metric_type_1 == 3)
    {
        // BL
        if (metric_type_2 == 1)
        {
            double XKS[NDIM];
            X_MKSToKS_v2(X1, XKS);
            X_KSToBL(XKS, X2);
            return;
        }

        // KS
        else if (metric_type_2 == 2)
        {
            X_MKSToKS_v2(X1, X2);
            return;
        }

        else
        {
            throw std::invalid_argument("Only MKS to BL and KS supported");
            return;
        }
    }

    else
    {
        throw std::invalid_argument("Invalid metric type arguments");
        return;
    }
}

void TransformFourVectors(int metric_type_1, int metric_type_2, double X1[NDIM], double T1[NDIM], double T2[NDIM])
{
    // BL to KS / MKS
    if (metric_type_1 == 1)
    {
        // KS
        if (metric_type_2 == 2)
        {
            T_BLToKS(T1, T2, X1);
            return;
        }

        // MKS
        else if (metric_type_2 == 3)
        {
            double XKS[NDIM], TKS[NDIM];
            TransformCoordinates(1, 2, X1, XKS);
            T_BLToKS(T1, TKS, X1);
            T_KSToMKS_v2(TKS, T2, XKS);
            return;
        }

        else
        {
            throw std::invalid_argument("Only BL to KS and MKS supported");
            return;
        }
    }

    // KS to BL / MKS
    else if (metric_type_1 == 2)
    {
        // BL
        if (metric_type_2 == 1)
        {
            T_KSToBL(T1, T2, X1);
            return;
        }

        // MKS
        else if (metric_type_2 == 3)
        {
            T_KSToMKS_v2(T1, T2, X1);
            return;
        }

        else
        {
            throw std::invalid_argument("Only KS to BL and MKS supported");
            return;
        }
    }

    // MKS to BL / KS
    else if (metric_type_1 == 3)
    {
        // BL
        if (metric_type_2 == 1)
        {
            double XKS[NDIM], TKS[NDIM];
            TransformCoordinates(3, 2, X1, XKS);
            T_MKSToKS_v2(T1, TKS, X1);
            T_KSToBL(TKS, T2, XKS);
            return;
        }

        // KS
        else if (metric_type_2 == 2)
        {
            T_MKSToKS_v2(T1, T2, X1);
            return;
        }

        else
        {
            throw std::invalid_argument("Only MKS to BL and KS supported");
            return;
        }
    }

    else
    {
        throw std::invalid_argument("Invalid metric type arguments");
        return;
    }
}

void UpperToLower(int metric_type, double X[NDIM], double Ucon[NDIM], double Ucov[NDIM])
{    
    double gcov[NDIM][NDIM];
    GcovFunc(metric_type, X, gcov);
    Ucov[0] = gcov[0][0]*Ucon[0] + gcov[0][1]*Ucon[1] + gcov[0][2]*Ucon[2] + gcov[0][3]*Ucon[3];
    Ucov[1] = gcov[0][1]*Ucon[0] + gcov[1][1]*Ucon[1] + gcov[1][2]*Ucon[2] + gcov[1][3]*Ucon[3];
    Ucov[2] = gcov[0][2]*Ucon[0] + gcov[2][1]*Ucon[1] + gcov[2][2]*Ucon[2] + gcov[2][3]*Ucon[3];
    Ucov[3] = gcov[0][3]*Ucon[0] + gcov[3][1]*Ucon[1] + gcov[3][2]*Ucon[2] + gcov[3][3]*Ucon[3];
    return;
}

// For the three matrices
void UpperToLower3(int metric_type, double X[NDIM], double Ucon3[NDIM-1], double Ucov3[NDIM-1])
{    
    double gcov[NDIM][NDIM];
    GcovFunc(metric_type, X, gcov);
    // mul;tiplying with gamma_ij
    Ucov3[0] = gcov[1][1]*Ucon3[0] + gcov[1][2]*Ucon3[1] + gcov[1][3]*Ucon3[2];
    Ucov3[1] = gcov[2][1]*Ucon3[0] + gcov[2][2]*Ucon3[1] + gcov[2][3]*Ucon3[2];
    Ucov3[2] = gcov[3][1]*Ucon3[0] + gcov[3][2]*Ucon3[1] + gcov[3][3]*Ucon3[2];
    return;
}


// Solve quadratic equation for u^0 given the spatial components of the 4 velocity and the ideal MHD b field
bool Get4Velocity(int metric_type, double X[NDIM], double ui_con[NDIM-1], double Ui_con[NDIM], double b_con[NDIM])
{

    bool success = false;
    double b_cov[NDIM];
    UpperToLower(metric_type, X, b_con, b_cov);
    int i,j;

    double gcov[NDIM][NDIM];
    GcovFunc(metric_type, X, gcov);
    double A=0., B=0., C=0.;
    A = gcov[0][0];
    for (i=0;i<NDIM-1;i++)
    {
        B += 2 * gcov[0][i+1] * ui_con[i];
        for (j=0; j<NDIM-1;j++)
        {
            C += gcov[i+1][j+1] * ui_con[i] * ui_con[j];
        }
        Ui_con[i+1] = ui_con[i];
    }
    C +=1.;
    // Obtain two roots for the answer
    double root1, root2;
    root1 = (-B + sqrt(SQR(B) - 4*A*C))/(2*A);
    root2 = (-B - sqrt(SQR(B) - 4*A*C))/(2*A);

    // Pick the root that gives b^mu u_mu = 0
    Ui_con[0] = root1;
    double norm = 0.;
    for (i=0; i<NDIM;i++)
    {
        norm += Ui_con[i] * b_cov[i];
    }

    // If the first root fails, try second one
    if (fabs(norm) > 1.e-2)
    {
        Ui_con[0] = root2;
        norm = 0.;
        for (i=0; i<NDIM;i++)
        {
            norm += Ui_con[i] * b_cov[i];
        }
        if (fabs(norm) > 1.e-2)
        {
            success = false;
            print("WARNING: Roots don't satisfy u^mu b_mu = 0");
            return success;
        }
        else
        {
            success = true;
            return success;
        }
    }

    else
    {
        success = true;
        return success;
    }
}

// Tests for 4 velocity satisfying norm, magnetic field and lorentz factor test
void Test4Velocity( int id, int metric_type, double X[NDIM], double U_con[NDIM], 
                    double b_con[NDIM], double _lfac, double alpha)
{
    // Norm test
    double norm_vel = 0.;
    norm_vel = Norm_con(metric_type, X, U_con);

    // b^mu test
    double b_cov[NDIM], norm_bmu=0.;
    UpperToLower(metric_type, X, b_con, b_cov);
    int i;
    for (i=0;i<NDIM;i++)
    {
        norm_bmu += U_con[i] * b_cov[i];
    }

    // lfac test
    double lfac = U_con[0] * alpha;

    print("******** Begin ***************");
    printvar("Test ", (double)id);
    printvar("Norm of 4 vel", norm_vel);
    printvar("Dot product with bmu", norm_bmu);
    printvar("lfac test: percent diff", 100*fabs(1 - (lfac/_lfac)));
    print("******** End *****************");
    return;
}

// Helper functions to make a tetrad
void normalize(double *vcon, double Gcov[NDIM][NDIM])
{
	int k, l;
	double norm;

	norm = 0.;
	for (k = 0; k < 4; k++)
		for (l = 0; l < 4; l++)
			norm += vcon[k] * vcon[l] * Gcov[k][l];

	norm = sqrt(fabs(norm));
	for (k = 0; k < 4; k++)
		vcon[k] /= norm;

	return;
}


void project_out(double *vcona, double *vconb, double Gcov[NDIM][NDIM])
{

	double adotb, vconb_sq;
	int k, l;

	vconb_sq = 0.;
	for (k = 0; k < 4; k++)
		for (l = 0; l < 4; l++)
			vconb_sq += vconb[k] * vconb[l] * Gcov[k][l];

	adotb = 0.;
	for (k = 0; k < 4; k++)
		for (l = 0; l < 4; l++)
			adotb += vcona[k] * vconb[l] * Gcov[k][l];

	for (k = 0; k < 4; k++)
		vcona[k] -= vconb[k] * adotb / vconb_sq;

	return;
}

void normalize_null(double Gcov[NDIM][NDIM], double K[])
{
	int k, l;
	double A, B, C;

	/* pop K back onto the light cone */
	A = Gcov[0][0];
	B = 0.;
	for (k = 1; k < 4; k++)
		B += 2. * Gcov[k][0] * K[k];
	C = 0.;
	for (k = 1; k < 4; k++)
		for (l = 1; l < 4; l++)
			C += Gcov[k][l] * K[k] * K[l];

	K[0] = (-B - sqrt(fabs(B * B - 4. * A * C))) / (2. * A);

	return;
}

double delta(int i, int j)
{
	if (i == j)
		return (1.);
	else
		return (0.);
}

void lower(double *ucon, double Gcov[NDIM][NDIM], double *ucov)
{

	ucov[0] = Gcov[0][0] * ucon[0]
	    + Gcov[0][1] * ucon[1]
	    + Gcov[0][2] * ucon[2]
	    + Gcov[0][3] * ucon[3];
	ucov[1] = Gcov[1][0] * ucon[0]
	    + Gcov[1][1] * ucon[1]
	    + Gcov[1][2] * ucon[2]
	    + Gcov[1][3] * ucon[3];
	ucov[2] = Gcov[2][0] * ucon[0]
	    + Gcov[2][1] * ucon[1]
	    + Gcov[2][2] * ucon[2]
	    + Gcov[2][3] * ucon[3];
	ucov[3] = Gcov[3][0] * ucon[0]
	    + Gcov[3][1] * ucon[1]
	    + Gcov[3][2] * ucon[2]
	    + Gcov[3][3] * ucon[3];
	return;
}


/* make orthonormal basis 
   first basis vector || U
   second basis vector || B
*/
void make_tetrad(double Ucon[NDIM],
		 double Gcov[NDIM][NDIM], double Econ[NDIM][NDIM],
		 double Ecov[NDIM][NDIM])
{
	int k, l;
	double norm;

	/* econ/ecov index explanation:
	   Econ[k][l]
	   k: index attached to tetrad basis
	   index down
	   l: index attached to coordinate basis 
	   index up
	   Ecov[k][l]
	   k: index attached to tetrad basis
	   index up
	   l: index attached to coordinate basis 
	   index down
	 */

	/* start w/ time component parallel to U */
	for (k = 0; k < 4; k++)
		Econ[0][k] = Ucon[k];
	normalize(Econ[0], Gcov);

	/*** done w/ basis vector 0 ***/

    // Trial vector is along the rhat direction
    double trial[NDIM] = {0., 1., 0., 0.};

	/* now use the trial vector in basis vector 1 */
	/* cast a suspicious eye on the trial vector... */
	norm = 0.;
	for (k = 0; k < 4; k++)
		for (l = 0; l < 4; l++)
			norm += trial[k] * trial[l] * Gcov[k][l];
	if (norm <= SMALL_VECTOR) {	/* bad trial vector; default to radial direction */
		for (k = 0; k < 4; k++)	/* trial vector */
			trial[k] = delta(k, 1);
	}

	for (k = 0; k < 4; k++)	/* trial vector */
		Econ[1][k] = trial[k];

	/* project out econ0 */
	project_out(Econ[1], Econ[0], Gcov);
	normalize(Econ[1], Gcov);

	/*** done w/ basis vector 1 ***/

	/* repeat for x2 unit basis vector */
	for (k = 0; k < 4; k++)	/* trial vector */
		Econ[2][k] = delta(k, 2);
	/* project out econ[0-1] */
	project_out(Econ[2], Econ[0], Gcov);
	project_out(Econ[2], Econ[1], Gcov);
	normalize(Econ[2], Gcov);

	/*** done w/ basis vector 2 ***/

	/* and repeat for x3 unit basis vector */
	for (k = 0; k < 4; k++)	/* trial vector */
		Econ[3][k] = delta(k, 3);
	/* project out econ[0-2] */
	project_out(Econ[3], Econ[0], Gcov);
	project_out(Econ[3], Econ[1], Gcov);
	project_out(Econ[3], Econ[2], Gcov);
	normalize(Econ[3], Gcov);

	/*** done w/ basis vector 3 ***/

	/* now make covariant version */
	for (k = 0; k < 4; k++) 
    {
		/* lower coordinate basis index */
        lower(Econ[k], Gcov, Ecov[k]);
	}

	/* then raise tetrad basis index */
	for (l = 0; l < 4; l++) 
    {
		Ecov[0][l] *= -1.;
	}
}

/* input and vectors are contravariant (index up) */
void coordinate_to_tetrad(double Ecov[NDIM][NDIM], double K[NDIM],
			  double K_tetrad[NDIM])
{
	int k;

	for (k = 0; k < 4; k++) {
		K_tetrad[k] =
		    Ecov[k][0] * K[0] +
		    Ecov[k][1] * K[1] +
		    Ecov[k][2] * K[2] + Ecov[k][3] * K[3];
	}
}

/* input and vectors are contravariant (index up) */
void tetrad_to_coordinate(double Econ[NDIM][NDIM], double K_tetrad[NDIM],
			  double K[NDIM])
{
	int l;

	for (l = 0; l < 4; l++) {
		K[l] = Econ[0][l] * K_tetrad[0] +
		    Econ[1][l] * K_tetrad[1] +
		    Econ[2][l] * K_tetrad[2] + Econ[3][l] * K_tetrad[3];
	}

	return;
}


// Get the norm starting with a covariant vector
double Norm_cov(int metric_type, double X[NDIM], double V_cov[NDIM])
{
	double norm=0.;
	double gcon[NDIM][NDIM];
	int i,j;
	GconFunc(metric_type, X, gcon);
	for (i=0;i<NDIM;i++)
	{
		for (j=0;j<NDIM;j++)
		{
			norm += gcon[i][j] * V_cov[j] * V_cov[i];
		}
	}
	return norm;
}

// Get the norm starting with a contravariant vector
double Norm_con(int metric_type, double X[NDIM], double V_con[NDIM])
{
	double norm=0.;
	double gcov[NDIM][NDIM];
	int i,j;
	GcovFunc(metric_type, X, gcov);
	for (i=0;i<NDIM;i++)
	{
		for (j=0;j<NDIM;j++)
		{
			norm += gcov[i][j] * V_con[j] * V_con[i];
		}
	}
	return norm;
}

// Dot product of two contravariant vectors
// vec_type 1: both are contravariant, 2: both are covariant
double ADotB(int metric_type, int vec_type, double X[NDIM], double a[NDIM], double b[NDIM])
{
    double dot = 0.;
    double g[NDIM][NDIM];
    if (vec_type == 1)
    {
        GcovFunc(metric_type, X, g);
    }

    else if (vec_type == 2)
    {
        GconFunc(metric_type, X, g);
    }

    else
    {
        std::cout<<"vec type must be 1 or 2"<<std::endl;
        return dot;
    }

    int i,j;
    for (i=0; i<NDIM; i++)
    {
        for (j=0; j<NDIM; j++)
        {
            dot += g[i][j] * a[i] * b[j];
        }
    }   
    return dot;
}

// Get the determinant of the covariant metric
double Detgcov(int metric_type, double X[NDIM])
{
    double det = 0.;
    double g[NDIM][NDIM];
    GcovFunc(metric_type, X, g);

    det += (g[0][0]*g[1][1]*g[2][2]*g[3][3] + g[0][0]*g[1][2]*g[2][3]*g[3][1] + 
            g[0][0]*g[1][3]*g[2][1]*g[3][2] - g[0][0]*g[1][3]*g[2][2]*g[3][1]);

    det += (-g[0][0]*g[1][1]*g[2][3]*g[3][2] - g[0][1]*g[1][0]*g[2][2]*g[3][3] - 
            g[0][2]*g[1][0]*g[2][3]*g[3][1] - g[0][3]*g[1][0]*g[2][1]*g[3][2]);

    det += (g[0][3]*g[1][0]*g[2][2]*g[3][1] + g[0][2]*g[1][0]*g[2][1]*g[3][3] + 
            g[0][1]*g[1][2]*g[2][0]*g[3][3] + g[0][2]*g[1][3]*g[2][0]*g[3][1]);

    det += (g[0][3]*g[1][1]*g[2][0]*g[3][2] - g[0][3]*g[1][2]*g[2][0]*g[3][1] - 
            g[0][2]*g[1][1]*g[2][0]*g[3][3] - g[0][1]*g[1][3]*g[2][0]*g[3][2]);

    det += (-g[0][1]*g[1][2]*g[2][3]*g[3][0] - g[0][2]*g[1][3]*g[2][1]*g[3][0] - 
            g[0][3]*g[1][1]*g[2][2]*g[3][0] + g[0][3]*g[1][2]*g[2][1]*g[3][0]);

    det += (g[0][2]*g[1][1]*g[2][3]*g[3][0] + g[0][1]*g[1][3]*g[2][2]*g[3][0]);

    return det;
}


// Return the determinant of the 3 metric
double Detgammacov(int metric_type, double X[NDIM])
{
    double gcon[NDIM][NDIM];
    GconFunc(metric_type, X, gcon);
    double alpha = sqrt(-1/gcon[0][0]); 

    double detg = Detgcov(metric_type, X);
    double detgamma = -detg/SQR(alpha);
    return detgamma;
}

// Returns the sign of a variable
double sgn(int x)
{
    if (x > 0)
    {
        return 1.;
    }
    else if (x < 0)
    {
        return -1.;
    }
    else
    {
        return 0.;
    }
}

// Compute Levi-Civita permutation without metric given indices
double GetEta(int a1, int a2, int a3, int a4)
{
    double lc = sgn(a2-a1)*sgn(a3-a2)*sgn(a3-a1)*sgn(a4-a3)*sgn(a4-a2)*sgn(a4-a1);
    return lc;
}

double GetEta(int a1, int a2, int a3)
{
    double lc = sgn(a2-a1)*sgn(a3-a2)*sgn(a3-a1);
    return lc;
}
