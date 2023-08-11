struct prop // time-dependent stuff...
{
	alfa, dt, kmax, kout;
}

struct params
{
	di, wi;
}

struct grid
{
	xmin, xmax, dx, am, pi; // -------- - autocorrelation function and probability;
}


struct prop prop;
struct params params;
struct grid grid;
//struct grid grid;

void plot(iwr, number_of_points, psi, cut);
void psigaus(number_of_points, psi, psib, v, alfa, xi, pxi, q0, nb); // gaussian localized wavefunction, may be moving (exponential factor)
void correl(number_of_points, cr, ham, eig, q0, nb); //supposedly correlation function from Eq 21 of Garashchuk
void aver(number_of_points, psi0, psi, c); // integrates a product of psi and complex conjugate of psi0, puts result in c - Griffith's c_n?
void ham_dvr(number_of_points, v, psi, cr, ham, w) // Miller-Colbert, why send psi here if it is not used inside;
void pes0(np, x, v); // creating potential v, currently symmetric fourth-degree polynomial

//all stuff above this point was not there in the fortran file. I guess converter created it

program qm()
{
	//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
	//c--------------------- MODEL LANDSCAPE --------------------------------
	//cccccccc use Gaussians instead of the eigenstates ccccccccccccccccccccc
	//ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

	// implicit double(a - h, o - z)

	int*4 number_of_points; // number of points

	const nx = 2048;
	const nbx = 50;
	double v(nx), iso, eig(nx), ham(nx, nx), q0(0:nbx); // this syntax with : means beginning and end. Probably to force 0-index
	//arrys in Fortran start from one by default . In c it does not matter much and syntax does not make sense.
	// code converter is just not smart enough
	double _Complex im, c, psi(nx), psib(nx), cr(nx);

	//COMMON/prop/alfa, dt, kmax, kout
	//COMMON/params/di, wi
	//COMMON/grid/xmin, xmax, dx, am, pi
	//c---------- some constants  ---------------------
	im = (0.d0, 1.d0); // imaginary one
	pi = 4.d0*atan(1.d0); // just pi, for some reason it is not just using numerical value

	// openumber_of_pointsng a whole bunch of files...
	fopen(15, file = "proj");
	fopen(16, file = "qwp");
	fopen(17, file = "cor");
	fopen(18, file = "eig") //eigenvalues;
	fopen(20, file = "wf0");
	fopen(21, file = "wft");
	fopen(22, file = "wfg");
	fopen(54, file = "vec"); //eigenvectors;
	fopen(77, file = "pop");
	fopen(15, file = "IN", status = "old"); // INPUT FILE with some parameters
	//c time propagation parameters - "15"; refers to file called IN

	scanf(15, *) dt, kmax, kout; // reads line from IN file, current values are: 50d0, 1000, 50 respectively
									// wave packet parameters
	scanf(15, *) alfa, starting_x_far, momentum; // 0.65d0, -17.0d0, 0.0d0;
								   //localized wavefunction parameters
	scanf(15, *) xmin; //lower left grid - 6.0d0;
	scanf(15, *) am; // mass = 1;
	scanf(15, *) number_of_points; //512 number of steps in grid;

	// guess just printing parametetrs to scree to let user know what they were.
		printf(*, *) "MASS", am;
	printf(*, *);

	printf(*, *) "GWP", alfa, starting_x_far, momentum;
	printf(*, *);

	xmax = -xmin; // from - 6 to 6...;
	dx = abs(xmax - xmin) / (number_of_points - 1) //step of the wavefunction argument;

	printf(*, *) "GRID:", dx, xmin, xmax;
	printf(*, *); -------------------------------------------------- - ;
	psigaus(number_of_points, psi, psib, v, alfa, starting_x_far, momentum, q0, nb); 
	// create gaussian wavefunction; starting_x_far = -17 for some reason; the two wells are at +/- 2.32

	ham_dvr(number_of_points, v, psi, cr, ham, eig) // do matrix diagonalization;

	//looks like main program is unfinished, as it does not go beyond Miller-Colbert (no decomposition of localized states)
	//need to check if this last subroutine calls on something else from inside... Does not look like that...

	//g1000:
	// format(20(e14.7, 1x))

	exit(0);
} //----end of main function - ---------------------------------------------------------------------------;

// --------------------------------------------
void plot(iwr, number_of_points, psi, cut)
{
	// implicit double(a - h, o - z)  - this is something that converter just did not know how to convert.
	//maybe necessary, maybe not
	//c----- write down wavefunction --------------
	const nx = 2048;
	double _Complex psi[nx];
	//COMMON/grid/xmin, xmax, dx, am, pi
	for (i = 0; i < number_of_points; i++) // change indexes to scheme starting from zero
	{
		x = grid.xmin + grid.dx*i;
		printf(iwr, *) x, abs(psi[i]) // looks like it plots the wavefunction;
	}
	return;
}

//c------------------Gaussian wavefunction-----------
void psigaus(number_of_points, psi, psib, v, alfa, xi, pxi, q0, nb)
{
	// alfa is inverse of STD squared. = 0.65
	// number_of_points is number of points
	// psi, psib - Gaussian wavepacket wavefunction and its complex conjugate, caluclated inside
	// v - potential, used only for some checking purposes
	// xi - center of the Gaussian, currently starting_x_far outside - currently -17
	//pxi - momentum for free-particle part, currently 0
	//q0 - defined below as an array of size nb, does not seem to be used
	//nb - not defined anywhere... number of barriers?
	const nx = 2048;
	const nbx = 50;
	double _Complex im, psi[nx], psib[nx];
	double v[nx], work[50], x[nx], q0[0:nbx];

	//COMMON/grid/xmin, xmax, dx, am, pi

	im = (0d0, 1d0) // imaginary one, i(not to be confused with cycle variable i);
		an = dsqrt(dsqrt(alfa / grid.pi)) // amplitude of the wavefunction;

	for (i = 0; i < number_of_points; i++) // to 512 ; define x-grid
	{
		x[i] = grid.xmin + grid.dx*i;
	} //end for

	pes0(number_of_points, x, v) // define potential v[i] = aq*x[i]**4 + bq*x[i]**2 + cq // double well turns array v with 512 elements;

	pot = 0d0;
	anrm = 0d0;

	for (i = 0; i < number_of_points, i++) //define WF; START LOOP
	{
		qx = x[i] - xi //displacement from the far point (-17)
		//xi = starting_x_far = -17 - why, given that minima of the double - well potential are at + / - 2.32;
		// // most likely this is assuming that we do not have outside walls...
			//alfa = 0.65
		psi[i] = an*exp(-alfa*qx**2 / 2 + im*pxi*qx) // pxi = 0 so at the moment the second term is turned off
			// so the wavefunction is currently real.
			//that second term looks like free-particle wavefunction exp(i p_x x)
		psib[i] = conjg(psi[i]); //complex conjugate of the wavefunction
		// if second term is on, complex conjugate has same exponential factor with opposite sign and it cancells out

		pot = pot + abs(psi[i])**2 * v[i] 
			// looks like integral of wavefunction squared times potential. Not sure why, does not do anything with it later
		anrm = anrm + abs(psi[i])**2 
			// looks like sum of wavefunction squared, not taking into accoung step size; the latter is done before print
		rho = abs(psi[i])**2 * grid.dx; // looks like current wavefunction squared, times the step size.;

		if (rho < e-16) then rho = e-16;

		printf(20, *) x[i], dsqrt(rho), v[i]; // printing results into file 20, amplitude of the wavefunction (see above)

	} // end for
	// created a Gaussian non-stationary wavefunction, possibly with movement
	// 
	// the following is some diagnostics, I guess
	printf(*, *) "nrm", anrm*grid.dx;
	printf(*, *) "pot", pot*grid.dx; //OK, it does print that "pot", but does not use it in any calculations...
	printf(*, *) "kin", alfa / grid.am / 2d0;

	return;
} //end of subroutine generating the localized MOVING wavefunction


//------------------------------------------------------- - -----------------------;
// --------------------------------------------
void correl(int number_of_points, double * cr, double * ham, double * eig, double q0, int nb)
//time-dependent wavefunction psi is re-generated inside from eigenvectors and coefficients cr(?)
// q0 looks like the list of barrier walls
//nb - number of barriers?
{
	// implicit double(a - h, o - z)
	const nx = 2048;
	const nbx = 50;
	double _Complex psi[nx], ot[0:10], c, tc, im, ct[nx], z; // cr[nx]
	double yr[3], ov[0:10], f[0:10]; // eig[nx], ham[nx, nx]  // again, 0:XXX is uncorrected Fortran syntax, meaning that indexes start from 0
	double q0[0:nbx], pop[0:nbx], x[nx];
	//COMMON/prop/alfa, dt, kmax, kout
	//COMMON/grid/xmin, xmax, dx, am, pi - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - autocorrelation function and probability

	printf(*, *) "COMPUTING COR FUN", prop.dt, prop.kmax, prop.kout;

	im = (0d0, 1d0); //imaginary one

	for (j = 0; j < number_of_points; j++) // define x-points
	{
		x[j] = grid.xmin + grid.dx*j;
	}

	q0[nb] = grid.xmax //define RHS grid for binnumber_of_pointsng;
		// if nb is number of barriers, q0 last barrier is located at the end of the grid

	for (k = 0; k <= prop.kmax; k++) // looks like cycle over time; kmax=1000 
	{
		t = k*prop.dt; // dt is currently 50 of unknown time units
		tc = -im*t; // shortcut for exp (-iE t/h_bar), h_bar=1 here
		kz = mod(k, prop.kout); //kout tells how often to produce output
		if (kz == 0) printf(*, *) "WF:", t, k, prop.kout;

		c = (0d0, 0d0);

		for (i = 0; i<nb; i++)
		{
			pop[i] = 0d0;
		} // end of for

		for (j = 0; j< number_of_points; j++)
		{
			c = c + abs(cr[j])**2 * exp(eig[j] * tc); // not sure what this is... not doing anything with it anyway
			ct[j] = cr[j] * exp(tc*eig[j]); // time dependence is attached to the coefficients cr that came from outside
		}
			//compute the wavefunction at t - ----------------;
		wf_nrm = 0d0;
		for (j = 0; j< number_of_points; j++) // number of points 512
		{
			psi[j] = (0d0, 0d0); // wavefunction reconstructed as linear superposition of eigenvecotrs

			for (i = 0; i < number_of_points; i++) // computed time - dependent WF
			{
				psi[j] = psi[j] + ct[i] * ham[j, i]; // time-dependent localized wavefunction
			}

			rho = abs(psi[j])**2;
			wf_nrm = wf_nrm + rho // normalization;
			for (i < nb; i = 0; i--) // populations, cycle over wells/barirers
			{
				if (x[j] < q0[i + 1]) && (x[j] > q0[i]) pop[i] = pop[i] + rho; 
				// populations rho_n, how much of the wavefunction fits into well n
				// if q0[nb] is the right wall at the end of the grid, then the above range of x corresponds to contents of the last well.
			} //end of i cycle

			if (kz == 0)
			{
				z = psi[j] / dsqrt(grid.dx);
				printf(21, 1000) x[j], z;  // format(20(e14.7, 1x));
			}

		}// end of j cycle

		printf(77, *) t, (pop[i], i = 0, nb); // how much given state contributes to the well

	} // end of time cycle

	return;
} //end of subroutine; so far it just produces time-dependent wavefunction 

//beginning of some averaging subroutine...;

// --------------------------------------------
void aver(number_of_points, double _Complex * psi0, double _Complex * psi, double _Complex c) 
// looks like it is multiplying psi with complex conjugate of psi0 and sums up.
//numerical integration. Calculates correlation function between moving and static wavefunctions
{
	// implicit double(a - h, o - z)
	const nx = 2048;

	c = (0.d0, 0.d0);
	for (i = 0; i < number_of_points; i++)
	{
		c = c + psi[i] * conjg(psi0[i]);
	}
	c = c*grid.dx; // i guess multiplying by x-interval is necessary for proper normalization of the result
	return; // returns the value of correlation function for one particular time (one particular psi)
}


// end of averaging subroutine = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = = ;

// --------------------------------------------
void ham_dvr(number_of_points, v, psi, cr, ham, w) // why send psi here if it is not used inside
{
	//--set up Hamiltonian matrix using Miller - Colbert DVR - -------;

	const nx = 2048;
	const lwork = nx * 10;
	const eps = 1d - 6;
	double work[lwork], w[nx], ham[nx, nx], hsv[nx, nx], v[nx], ps0(nx), z(nx);
	double _Complex psi[nx], im, wfc[nx], cr[nx], zz;

	//COMMON/grid/xmin, xmax, dx, am, pi

//g1000:
	// format(81(e14.7, 1x))

	im = (0d0, 1d0);
	akx = 0.5d0 / grid.am / grid.dx / grid.dx;
	p3 = grid.pi*grid.pi / 3d0; //Miller - Colbert DVR;

	n1 = number_of_points; //512;
			//c----------------- rewrite potential as a vector -----------
			//c----------------- construct the hamiltonumber_of_pointsan matrix --------------
	for (j = 0; j < number_of_points; j++) // to 512)\n
		{
			for (i = 0; i < number_of_points; i++)
			{
				ham[i, j] = 0d0;
				if (i == j) ham[i, j] = v[i] + akx*p3 // includes both kinetic and potential energy;
					if (i != j)
					{
						ham[i, j] = akx*(-1d0)**(j - i) * 2d0 / (j - i)**2;
					} //end of if
				hsv[i, j] = ham[i, j];
			} // end of for
		} //end of for

	printf(*, *) "DIAGONALIZATION";

	dsyev("V", "L", n1, ham, nx, w, work, lwork, info) // finding eigenvalues and eigenvectors;
													   //DSYEV(JOBZ, UPLO, N, A, LDA, W, WORK, LWORK, INFO)
													   // V save eigenvalues and eigenvectors
													   // lower triangle is stored
													   // N = n1 the order of the matrix
													   // ham
													   //(input/output) double array, //DIMENSION(LDA, N)
													   //*          On entry, the symmetric matrix A.  If UPLO = 'U', the
													   //*          leading N-by-N upper triangular part of A contains the
													   //*          upper triangular part of the matrix A.  If UPLO = 'L',
													   //*          the leading N-by-N lower triangular part of A contains
													   //*          the lower triangular part of the matrix A.
													   //*          On exit, if JOBZ = 'V', then if INFO = 0, A contains the
													   //*          orthonormal eigenvectors of the matrix A.
													   //*          If JOBZ = 'N', then on exit the lower triangle (if UPLO='L')
													   //*          or the upper triangle (if UPLO='U') of A, including the
													   //*          diagonal, is destroyed.eading dimension of the matrix INFO = 0, the eigenvalues in ascending order.

		printf(*, *) "INFO", info, work[1];

	printf(*, *) "EIGENVALS", w[1], w[2], w[3];

	for  i = 0 to  49 //print first 50 eigenvalues
	{
		printf(18, *) i, w[i] //write eigenvalues to file;
	}


	//c---------------------- check the eigenvector -----------------
	printf(*, *) "ENTER VECTOR TO TEST";
	// at this point ham contains orthonormal eigfenvectors

	nt = 2;
	eig = w[nt]; // eigenvalue returned by dsyev;
	dev = 0d0;
	dds = dsqrt(grid.dx);

	for (i = 0; i<= number_of_points-1, i++) // 30 is label)\n
	{
		ps0(i) = ham[i, nt] / dds;
		if (ps0(i) < eps) ps0(i) = eps;
		x = grid.xmin + i*grid.dx;
		printf(54, 1000) x, (ham[i, k] / dds, k = 1, 20); // write eigenvectors file; first column is position;

			z(i) = 0d0;

		for (j = 0; j< number_of_points; j++) // 31 is label)\n
		{
			z(i) = z(i) + hsv[i, j] * ham[j, nt] / eig; // hsv is backuped initial matrix
			//multiplying the original hamiltonian by second eigenvector and dividing by second eigenvalue
			// does it give back the second eigenvector?
		}
		dev = dev + abs(z(i) - ham[i, nt])**2; // deviation
			//------------------------------------------------------------------------TEST the eigenvector - ----------------------------------------------------------------------;
	} // end of for
	printf(*, *) "EIGENVECTOR TEST", dev;

		anv = 0d0;
		for (n = 0; n < n1; n++)
		{
			anv = anv + abs(ham[n, nt])**2; // normalization of one parricular eigenvector
		}

	printf(*, *) "VECTOR INDEX and NORM=", nt, anv;

	return;
}  //c--------------DVR thing ended-----------
//c--------------creating the potential--------------------
// --------------------------------------------
void pes0(int np, double * x, double * v) //c----- double well Batista model I pacs.jpca.5b12192 ---
{
		//c--------------------------------------------------------
	const aq = 0.04614589486d0;
	const bq = -0.5d0;
	const cq = 0d0;
	const nx = 2048;
	//double x[nx], v[nx]; // declaring arrays of maximal size and then not fully using them

	for (i = 0; i < np; i++) // to number of points = 512
	{
		v[i] = aq*x[i] * *4 + bq*x[i] * *2 + cq; // double well;
		//v[i] = x[i] * *2 / 2d0; // Harmonic Oscillator test;
	} //end of for
	return;
} //end of sub



