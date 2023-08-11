
// MitacsGUIDlg.cpp: 实现文件
//

#include "pch.h"
#define mainfile 1
#include "framework.h"
#include "MitacsGUI.h"
#include "MitacsGUIDlg.h"
#include "afxdialogex.h"
#include <iostream>   // std::cout
//#include <string>     // included in Glbal_variable.h
#include <time.h>
#include "Global_variables.h"
//#include <sstream>
#include <fstream>

using namespace std;


#ifdef _DEBUG
#define new DEBUG_NEW
#endif


// CAboutDlg dialog box for the "About" menu item of the application
class CAboutDlg : public CDialogEx
{
public:
	CAboutDlg();

// 对话框数据
#ifdef AFX_DESIGN_TIME
	enum { IDD = IDD_ABOUTBOX };
#endif

	protected:
	virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV data exchange and validation

// 实现
protected:
	DECLARE_MESSAGE_MAP()
};

CAboutDlg::CAboutDlg() : CDialogEx(IDD_ABOUTBOX)
{
}

void CAboutDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialogEx::DoDataExchange(pDX);
}

BEGIN_MESSAGE_MAP(CAboutDlg, CDialogEx)
END_MESSAGE_MAP()


// CMitacsGUIDlg 对话框


//some default values
CMitacsGUIDlg::CMitacsGUIDlg(CWnd* pParent /*=nullptr*/)
	: CDialogEx(IDD_MITACSGUI_DIALOG, pParent)
	, mass_text(_T("1.67e-26")) // I made it ten times the mass of the proton, to ensure that md^2 will have correct order, 10-46 kg m2
	, inhomog_text(_T("14700"))
	, inhFWHM_text(_T("70"))
	, beam_d_text(_T("0.5"))
	, crossection_text(_T("2.258e-16"))
	, homogRT_text(_T("300"))
	, homogfinalT_text(_T("0.5"))
	, S_text(_T("0.3"))
	, RT_text(_T("300"))
	, iniT_text(_T("300"))
	, finalT_text(_T("5"))
	, cool_t_text(_T("120"))
	, burn_freq_text(_T("14720"))
	, burnP_text(_T("1"))
	, burn_t_text(_T("60"))
	, scanN_text(_T("100"))
	, t_perscan_text(_T("10"))
	, min_scan_text(_T("14600"))
	, max_scan_text(_T("14800"))
	, scan_step_text(_T("1"))
	, recovery_t_text(_T("6000"))
	, recovery_points_text(_T("100"))
	, cycling_maxT_text(_T("60"))
	, Ncomplex_text(_T("1"))
	, Nwells_text(_T("10"))
	, curvature_text(_T("0"))
	, barrier_mean_text(_T("1200"))
	, barrier_STD_text(_T("50"))
	, lifetime_text(_T("3"))
	, stretch_text(_T("0.25"))
	, bottom_mean_text(_T("0"))
	, bottom_STD_text(_T("5"))
	, Nhires_text(_T("512"))
	, filename_text(_T("Not working yet"))
	, dist_text(_T("0.2")) // distance between minima of the wells, in nm. If rectangular barriers, this is one well 
	// plus one barrier, so barrier thickness is 0.1 nm or 1 A.
	, T_step_text(_T("5"))
	, BT_text(_T("5"))
	, HGK_step_text(_T("100"))
	, well_to_STD_ratio(_T("6"))
{
	m_hIcon = AfxGetApp()->LoadIcon(IDR_MAINFRAME);
}

void CMitacsGUIDlg::DoDataExchange(CDataExchange* pDX)
{
	CDialogEx::DoDataExchange(pDX);
	// writing data as strings into variables. 
	DDX_Text(pDX, IDC_EDIT1, mass_text);
	DDX_Text(pDX, IDC_EDIT2, inhomog_text);
	DDX_Text(pDX, IDC_EDIT3, inhFWHM_text);
	DDX_Text(pDX, IDC_EDIT4, beam_d_text);
	DDX_Text(pDX, IDC_EDIT5, crossection_text);
	DDX_Text(pDX, IDC_EDIT6, homogRT_text);
	DDX_Text(pDX, IDC_EDIT7, homogfinalT_text);
	DDX_Text(pDX, IDC_EDIT25, RT_text);
	DDX_Text(pDX, IDC_EDIT26, iniT_text);
	DDX_Text(pDX, IDC_EDIT27, finalT_text);
	DDX_Text(pDX, IDC_EDIT40, cool_t_text);
	DDX_Text(pDX, IDC_EDIT14, burn_freq_text);
	DDX_Text(pDX, IDC_EDIT15, burnP_text);
	DDX_Text(pDX, IDC_EDIT12, burn_t_text);
	DDX_Text(pDX, IDC_EDIT10, scanN_text);
	DDX_Text(pDX, IDC_EDIT11, t_perscan_text);
	DDX_Text(pDX, IDC_EDIT37, min_scan_text);
	DDX_Text(pDX, IDC_EDIT38, max_scan_text);
	DDX_Text(pDX, IDC_EDIT39, scan_step_text);
	DDX_Text(pDX, IDC_EDIT32, recovery_t_text);
	DDX_Text(pDX, IDC_EDIT42, recovery_points_text);
	DDX_Text(pDX, IDC_EDIT47, cycling_maxT_text);
	DDX_Text(pDX, IDC_EDIT17, Ncomplex_text);
	DDX_Text(pDX, IDC_EDIT18, Nwells_text);
	DDX_Text(pDX, IDC_EDIT22, curvature_text);
	DDX_Text(pDX, IDC_EDIT23, barrier_mean_text);
	DDX_Text(pDX, IDC_EDIT24, barrier_STD_text);
	DDX_Text(pDX, IDC_EDIT8, S_text);
	DDX_Text(pDX, IDC_EDIT41, lifetime_text);
	DDX_Text(pDX, IDC_EDIT33, stretch_text);
	DDX_Text(pDX, IDC_EDIT34, bottom_mean_text);
	DDX_Text(pDX, IDC_EDIT35, bottom_STD_text);
	DDX_Text(pDX, IDC_EDIT51, Nhires_text);
	DDX_Text(pDX, IDC_EDIT19, filename_text);
	DDX_Text(pDX, IDC_EDIT48, dist_text); // distance between minima, or barrier thickness plus well width
	DDX_Text(pDX, IDC_EDIT28, T_step_text); //this is called even before dialog appears, initial/default values of various parameters
	DDX_Text(pDX, IDC_EDIT50, BT_text);
	DDX_Text(pDX, IDC_EDIT13, HGK_step_text);
	DDX_Control(pDX, IDC_Equilibrium, Equilibrium);
	DDX_Control(pDX, IDC_Uniform, Uniform);
	DDX_Control(pDX, IDC_generate, generate);
	DDX_Control(pDX, IDC_balance, balance);
	DDX_Control(pDX, IDC_bottom, bottom);
	DDX_Control(pDX, IDC_barrier, barrier);
	DDX_Control(pDX, IDC_analytical, analytical);
	DDX_Control(pDX, IDC_scan, scan);
	DDX_Control(pDX, IDC_details, details);
	DDX_Control(pDX, IDC_one_step, one_step);
	DDX_Control(pDX, IDC_recovery_burn, recovery_burn);
	DDX_Text(pDX, IDC_EDIT16, well_to_STD_ratio);
}

BEGIN_MESSAGE_MAP(CMitacsGUIDlg, CDialogEx)
	ON_WM_SYSCOMMAND()
	ON_WM_PAINT()
	ON_WM_QUERYDRAGICON()
	//ON_LBN_SELCHANGE(IDC_LIST1, &CMitacsGUIDlg::OnLbnSelchangeList1)
	ON_BN_CLICKED(IDOK, &CMitacsGUIDlg::OnBnClickedOk)
	ON_BN_CLICKED(IDCANCEL, &CMitacsGUIDlg::OnBnClickedCancel)
	ON_BN_CLICKED(IDC_Equilibrium, &CMitacsGUIDlg::OnBnClickedEquilibrium)
	ON_BN_CLICKED(IDC_Uniform, &CMitacsGUIDlg::OnBnClickedUniform)
	ON_BN_CLICKED(IDC_generate, &CMitacsGUIDlg::OnBnClickedgenerate)
	ON_BN_CLICKED(IDC_balance, &CMitacsGUIDlg::OnBnClickedbalance)
	ON_BN_CLICKED(IDC_bottom, &CMitacsGUIDlg::OnBnClickedbottom)
	ON_BN_CLICKED(IDC_barrier, &CMitacsGUIDlg::OnBnClickedbarrier)
	ON_BN_CLICKED(IDC_analytical, &CMitacsGUIDlg::OnBnClickedanalytical)
	ON_BN_CLICKED(IDC_scan, &CMitacsGUIDlg::OnBnClickedscan)
	ON_BN_CLICKED(IDC_details, &CMitacsGUIDlg::OnBnClickeddetails)
	ON_BN_CLICKED(IDC_one_step, &CMitacsGUIDlg::OnBnClickedonestep)
	ON_BN_CLICKED(IDC_recovery_burn, &CMitacsGUIDlg::OnBnClickedrecoveryburn)
	ON_BN_CLICKED(IDC_Burn, &CMitacsGUIDlg::OnBnClickedBurn)
	ON_BN_CLICKED(IDC_SMSscan, &CMitacsGUIDlg::OnBnClickedSmsscan)
	ON_BN_CLICKED(IDC_recovery, &CMitacsGUIDlg::OnBnClickedrecovery)
	ON_BN_CLICKED(IDC_Thermocycling, &CMitacsGUIDlg::OnBnClickedThermocycling)
	//ON_EN_CHANGE(IDC_EDIT42, &CMitacsGUIDlg::OnEnChangeEdit42)
END_MESSAGE_MAP()


// CMitacsGUIDlg 消息处理程序

BOOL CMitacsGUIDlg::OnInitDialog()
{
	CDialogEx::OnInitDialog();

	// Add the "About..." menu item to the system menu.

		// IDM_ABOUTBOX must be within the system command range.


	ASSERT((IDM_ABOUTBOX & 0xFFF0) == IDM_ABOUTBOX);
	ASSERT(IDM_ABOUTBOX < 0xF000);

	CMenu* pSysMenu = GetSystemMenu(FALSE);
	if (pSysMenu != nullptr)
	{
		BOOL bNameValid;
		CString strAboutMenu;
		bNameValid = strAboutMenu.LoadString(IDS_ABOUTBOX);
		ASSERT(bNameValid);
		if (!strAboutMenu.IsEmpty())
		{
			pSysMenu->AppendMenu(MF_SEPARATOR);
			pSysMenu->AppendMenu(MF_STRING, IDM_ABOUTBOX, strAboutMenu);
		}
	}

	// Set the icon for this dialog box. When the main window of the application is not a dialog box, the frame will automatically do this


	SetIcon(m_hIcon, TRUE);			// 设置大图标
	SetIcon(m_hIcon, FALSE);		// 设置小图标

	// TODO: Add additional initialization code here

	

	return TRUE;  // Return unless the focus is set to the control TRUE
}

void CMitacsGUIDlg::OnSysCommand(UINT nID, LPARAM lParam)
{
	if ((nID & 0xFFF0) == IDM_ABOUTBOX)
	{
		CAboutDlg dlgAbout;
		dlgAbout.DoModal(); // this may require changing if stuff does not work
		//dlgAbout.Create(IDD_DIALOG1, this);
	}
	else
	{
		CDialogEx::OnSysCommand(nID, lParam);
	}
}

/// If you add a minimize button to the dialog box, you need the following code
// to draw the icon. For MFC applications that use the document/view model,
// This will be done automatically by the framework.


void CMitacsGUIDlg::OnPaint()
{
	if (IsIconic())
	{
		CPaintDC dc(this); // 用于绘制的设备上下文

		SendMessage(WM_ICONERASEBKGND, reinterpret_cast<WPARAM>(dc.GetSafeHdc()), 0);

		// 使图标在工作区矩形中居中
		int cxIcon = GetSystemMetrics(SM_CXICON);
		int cyIcon = GetSystemMetrics(SM_CYICON);
		CRect rect;
		GetClientRect(&rect);
		int x = (rect.Width() - cxIcon + 1) / 2;
		int y = (rect.Height() - cyIcon + 1) / 2;

		// 绘制图标
		dc.DrawIcon(x, y, m_hIcon);
	}
	else
	{
		CDialogEx::OnPaint();
	}
}

//当用户拖动最小化窗口时系统调用此函数取得光标
//显示。
HCURSOR CMitacsGUIDlg::OnQueryDragIcon()
{
	return static_cast<HCURSOR>(m_hIcon);
}



void CMitacsGUIDlg::OnLbnSelchangeList1()
{
	// TODO: 在此添加控件通知处理程序代码
}


void CMitacsGUIDlg::OnBnClickedCheck2()  // abort button
{
	// TODO: 在此添加控件通知处理程序代码
}

// initialization and cooling
void CMitacsGUIDlg::OnBnClickedOk() // generate energy landscapes and cool down button
{
	UpdateData(TRUE); //ensures that modified data from the interface gets into CString variables
	h_planck = 6.626068 * 1e-34;
	pi = 3.141592;
	energy_conversion_factor = 1e3;
	time_conversion_factor = 1e-3;
	mass = _ttof(mass_text);
	d_distance_of_two_wells = _ttof(dist_text); //in nm
	d_distance_of_two_wells = d_distance_of_two_wells * 1e-9; // in meters
	width_well = d_distance_of_two_wells / 2; //the width of a well or thickness of a barrier
	inhomogeneous_bandwidth = _ttof(inhFWHM_text);// cm^-1.  need conversion to MHz
	standard_deviation_of_translation = (inhomogeneous_bandwidth * 29.9792458 / 2.355) * energy_conversion_factor; // in MHz
	mean_value_of_translation = _ttof(inhomog_text); // in cm-1, need conversion to MHz
	mean_value_of_translation = mean_value_of_translation * 30000; //in MHz
	burn_frequency = _ttof(burn_freq_text);//frequency in cm-1, needs conversion to MHz.
	burn_frequency = burn_frequency * 30000; // in MHz
	burn_power = _ttof(burnP_text);//in mW. 1microW is 1e-3, as an example. 5, 10, and 15 W are taken from Kohler's paper.
	beam_diameter = _ttof(beam_d_text); //1/(2/sqrt(pi)); //cm. NB: In Kohler's paper the power is given in W/cm^2, so I have to conver the diameter in such a way that it gets rid of other constants in my equation (done on on September 12, 2012, in the second log book). In that calculation, area = pi r^2 = pi d^2/4. by assuming d = 2/sqrt(pi) the area would be 1cm^2. BTW sqrt(pi)/2 = 0.886 cm.
	absorption_cross_section_at_peak_RT = _ttof(crossection_text);//cm^2// September 12, 2012 (second log book)
	homogeneous_width_at_RT = _ttof(homogRT_text);//cm^-1 for Cp43
	homogeneous_width_at_BT = _ttof(homogfinalT_text);//cm^-1 for Cp43
	phonon_coupling = _ttof(S_text);// this is S or Huang–Rhys factor for Cp43.
	absorption_cross_section_at_peak_BT = absorption_cross_section_at_peak_RT * (homogeneous_width_at_RT / homogeneous_width_at_BT) * exp(-phonon_coupling);//=1.054e-12;//cm^2, see the notes on September 12, 2012 (in the second log book).
	photon_flux = burn_power * 1e-3 / (pi * (beam_diameter / 2) * (beam_diameter / 2) * h_planck * burn_frequency * 1e6); 
	//burn frequency in MHz by now, Planck constant in SI units

	time_interval_between_two_consecutive_burns = 1 / (absorption_cross_section_at_peak_BT * photon_flux); // in seconds
	//time_interval_between_two_consecutive_burns = time_interval_between_two_consecutive_burns * 1e6; // in mks
	whole_burning_time = _ttof(burn_t_text); //in seconds

	burn_time = _ttof(lifetime_text) * time_conversion_factor; //excited state lifetime in microseconds

	number_of_burn = (int)(whole_burning_time / time_interval_between_two_consecutive_burns); // number of acts of burn in HGK, both time intervals are in seconds
	HGK_measurement_frequency = _ttoi(HGK_step_text);//It means, how frequently the hole depth should be measured (e.g. every 100 acts of burn) 
	

	//cooling
	room_temperature = _ttof(RT_text);//Kelvin.// it could be the same as initial_temperature.
		//parameters related to the cooling for finding the distribution before the start of burning.
	final_temperature = _ttof(finalT_text); // Kelvin.
	initial_temperature = _ttof(iniT_text); // Kelvin.
	temperature_step = _ttof(T_step_text);// the step of reducing of temperature. NB. Changing of this parameter should not be done unless cooling time changes properly.
	whole_cooling_time = _ttof(cool_t_text); //in minutes
	whole_cooling_time = whole_cooling_time * 60 * 1e6;  // in microseconds
	cooling_time = whole_cooling_time * temperature_step / (room_temperature-final_temperature); // 59 steps
	//cooling time at each cooling step, in microseconds
	burn_temperature = _ttof(BT_text); // Kelvin.

	//Single Molecule Experiment
	number_of_scan = _ttoi(scanN_text); //one if we do not want to scan, only want the holeburn? Should we change it to zero?
	scan_time = _ttof(t_perscan_text);
	each_scan_time = scan_time * 1e6;// micro seconds. 
	minimum_frequency = _ttof(min_scan_text);
	minimum_frequency = minimum_frequency * 30000;
	maximum_frequency = _ttof(max_scan_text);
	maximum_frequency = maximum_frequency * 30000;
	frequency_step = _ttof(scan_step_text); // in cm-1
	frequency_step = frequency_step * 30000; // in MHz
	// in cm-1. Mehdi played with this value such that all wells are chosen once only.
	number_of_frequency_steps = (int)(((maximum_frequency - minimum_frequency) / frequency_step) + 1);
	time_interval_between_each_frequency_step = frequency_step * each_scan_time / (maximum_frequency - minimum_frequency); // in microseconds...

	//Landscapes

	ensemble_number = _ttoi(Ncomplex_text); // number of ensembles used. More ensembles are tested until we find 5000 with resonant wells.
	Hamiltonian_well_numbers = _ttoi(Nhires_text); // Should be 512/1024/2048 for optimal matrix solving - the size of the high resolution energy landscape used in Miller-Colbert algorithm mentioned in Garashchuk paper
	well_numbers = _ttoi(Nwells_text);
	point_numbers = well_numbers * 2 + 1;
	interval = width_well * point_numbers / 2; // as of July 22, 2021 this is determined by the well thickness
	//from center of the landscape to edge, in meters. Essentially it is a pseudo-shift of the parabola, allowing for loop 
	step = 2 * interval / point_numbers;
	mu = _ttof(barrier_mean_text);
	mu = mu * 30000; //in MHz
	sigma = _ttof(barrier_STD_text);
	sigma = sigma * 30000; //in MHz
	stretch = _ttof(stretch_text); //0.19; // ratio of excited and ground state BARRIERS
	// Bottoms of the wells in the ground state
	mean = _ttof(bottom_mean_text);
	mean = mean * 30000; //in MHz
	deviation = _ttof(bottom_STD_text);
	deviation = deviation * 30000; // in Mhz
	curvature_of_parabola = _ttof(curvature_text);
	well_STD_ratio = _ttof(well_to_STD_ratio);
	lambda_coefficient = width_well * sqrt(8 * mass * 1e6) * pi / h_planck;
	// lambda_coefficient is the constant part of lambda, which is calculated in the createRateMatrix.cpp; 

	//Recovery

	recovery_temperature = burn_temperature;//Kelvin.
	whole_recovery_time = _ttof(recovery_t_text); //in minutes
	whole_recovery_time = whole_recovery_time * 60 * 1000000; //in microseconds
	number_of_recovery_measurement = _ttoi(recovery_points_text);
	// it means that during the whole_recovery_time, we will stop the evolution of probabilities and calculate the spectum number_of_recovery_measurement (e.g. 60) times.
	recovery_time = whole_recovery_time / number_of_recovery_measurement;
	// it means that calculation for each step will take recovery_time.

	number_of_longtime_measurement = 1; // this is for simulations of ultra-long time evolution...
	whole_long_time = 13.0 * 365.0 * 24.0 * 60.0 * 60.0 * 1e9 * time_conversion_factor;
	// this time (13*365*24*60*60*1e9 = 13 years) is conidered as a longtime to calculate stationary state and then the recovery depth is compared with this distribution.
	long_time = whole_long_time / number_of_longtime_measurement;
	//const double long_time_step = long_time / 1.0;

	highest_experimental_cycling_temperature = _ttof(cycling_maxT_text);// This is the highest temperature during thermocycling experiment.
	highest_cycling_temperature = highest_experimental_cycling_temperature - burn_temperature;
	// Since the first temperature in thermocycling loop should be greater than burn temperature which results in having final temperature higher than what we have 
	//in practice (experimental_cycling_temperature).

	whole_thermocycling_time = whole_recovery_time;
	number_of_thermocycling_measurement = number_of_recovery_measurement;
	//const double cycle_temperature_step = cycle_temperature-burn_temperature / number_of_thermocycling_measurement;
	thermocycling_time = whole_thermocycling_time / number_of_thermocycling_measurement;// it means that calculation for each step will take thermocycling_time.

	all_frequency = well_numbers * ensemble_number;// number of all possible frequencies which come out of the whole set of molecules.

	total_bins = (int)sqrt(all_frequency); // number of bins for creating spectra. 
	// another option is to find minimal and maximal frequencies and divide it by frequency step
	averaged_probability_SMS_bin_number = 500;

	//Checkbox:
	generate_energy_landscape = 1;

	number_of_column = well_numbers * ensemble_number + 1; // for files containing probability records
	// the first column is for time so the number of column should be increased by 1
	// Code from constants.h ends
	// 
	// 
	// 
	//
	// Relevant code from ConsoleApplication4 (Bole Yi's latest version)
	// start and end are defined to calculate running time of program
	clock_t start, end; //long clock
	start = clock();
	//time_t now = time(0);
	//tm* localtm = localtime(&now);
	__time64_t longtime;
	errno_t err;
	struct tm  timeinfo;
	_time64(&longtime);
	err = _localtime64_s(&timeinfo, &longtime);
	//cout << "******* The program started at " << timeinfo.tm_hour << ":" << timeinfo.tm_min << ":" << timeinfo.tm_sec << " *************" << endl;
	srand((unsigned)time(0));
	// random number generator function in gsl library
	const gsl_rng_type* T;
	gsl_rng* r;
	gsl_rng_env_setup();

	T = gsl_rng_default;
	r = gsl_rng_alloc(T);
	gsl_rng_set(r, time(0)); //time(0) is used as seed for random number generator
	//
	//Allocating memory for landscape matrix and vector with excited-ground state baseline differences.
	//gsl_matrix* landscape_file = gsl_matrix_calloc(point_numbers, (3.0 * ensemble_number + 3));
	//Why 3-s ? At the first glance one needs only point_numbers times ensemble_number*2, for excited and ground states

	gsl_vector* translation = gsl_vector_calloc(ensemble_number); // shift between ground and excited state baselines (flat or parabola min)
	//one per molecule
	// not sure we care
	//
	int progression_percent;
	if (ensemble_number < 5)// for ensemble number smaller than 5, the progress indicator is pointless. 
	//After all calculation is done, there will be a message that 100% of calculation is done.
		progression_percent = 100;
	else// for bigger than 5 systems, the progression will be measured every 20%.
		progression_percent = 20;// 20%, when the calculation of 20% of ensemble_number is done, the program progression is 20%.
	int progression_value = (int)((progression_percent / 100.0) * ensemble_number);
	//
	//
	//

	starting_well = gsl_vector_calloc(ensemble_number); //resonant well for each molecule
	//
	int well_numbers_counter2 = 1;
	int produced_ensemble = 0; // counter of all ensembles produced, including those that do not have resonant wells
	// it shows how many systems has been produced so far total to end up with the ensemble_number system 
	//which have the same energy difference with respect to the burn energy.
	int binning=0;


	if (ensemble_number < 10)
		cutoff_energy_difference = 30 * energy_conversion_factor; // In MHz, i.e 100000 MHz or 1 cm-1
	else
		cutoff_energy_difference = 100 * energy_conversion_factor; // In MHz, 3 cm-1
	// it determines the difference between transition frequency and burn frequency for which we say that well is in resonance. 
	//In other words, this is almost as big as the bin step.
	int jj=0;
	gsl_vector* parabola = gsl_vector_calloc(point_numbers);
	for (double u = 0; u < 2 * interval; u += width_well)
	{
		gsl_vector_set(parabola, jj, curvature_of_parabola * (u - interval) * (u - interval));
		// + anharmonicity * (i-interval) * (i-interval) * (i-interval) );
		++jj;
	}

	energy_difference = gsl_matrix_calloc(well_numbers, ensemble_number);

	// for each well in each system the difference of ground and excited states is calculated
// as the difference of the stationary states with the largest PROJECTION rho for this well as of December 2022
// see GRASHCHUK'S paper about projections
// 
//
//variables that contain energy landscapes modulated by randomness
//
	excited_parabola_for_rate = gsl_matrix_calloc(ensemble_number, point_numbers);
	ground_parabola_for_rate = gsl_matrix_calloc(ensemble_number, point_numbers);
	// zig-zag landscapes with randomness; initially based on random numbers directly, 
	// eventually the bottoms are changed to most likely stationary states for the purpose of 
	// creating initial distributions A0_ground
	//
	eigenve = gsl_matrix_calloc(ensemble_number, Hamiltonian_well_numbers);
	//excited state eigenvalues - as of November 26, 2021 they are matrixes
	eigenvg = gsl_matrix_calloc(ensemble_number, Hamiltonian_well_numbers);
	//ground state eigenvalues
	// eigenvectors are not passed from generate... to other programs, they are used internally in generate_energy_landscape.
	//
	// 3D arrays converted to 1D arrays.
	//max array size 2GB, there are 6 such arrays, so 12 GB + several GB of smaller stuff. Good for 32GB RAM machine.
	//
	//We cannot have one variable that is a product of TE and attempt frequency. TE is a property of the barrier, and 
	// attempt frequency is property of the well, so the produc is different for R to L and L to R.
	//
	coeff_matrixe = gsl_vector_calloc(ensemble_number * well_numbers * Hamiltonian_well_numbers);
	//As of November 2022 these are coefficients rho describing projections of stationary states on particular wells. See Garashchuk's paper. NOT c_n of Griffiths.
	coeff_matrixg = gsl_vector_calloc(ensemble_number * well_numbers * Hamiltonian_well_numbers);
	//
	attempt_freq_matrix_e = gsl_vector_calloc(ensemble_number * well_numbers * Hamiltonian_well_numbers);
	attempt_freq_matrix_g = gsl_vector_calloc(ensemble_number * well_numbers * Hamiltonian_well_numbers);
	//
	//T-independent one-act tunneling probabilities T(E):
	TE_matrix_ground = gsl_vector_calloc(ensemble_number * well_numbers * Hamiltonian_well_numbers);
	TE_matrix_excited = gsl_vector_calloc(ensemble_number * well_numbers * Hamiltonian_well_numbers);
	//NB: no need for separate left to right and right to left matrices, one just needs to shift them by one, as these probabilities
	// actually correspond to barriers, not wells
	//
	v_barrier_height_excited = gsl_vector_calloc(ensemble_number * (well_numbers - 1));
	v_barrier_height_ground = gsl_vector_calloc(ensemble_number * (well_numbers - 1));
	// looks like these are used only for reporting barrier height distributions...	
	// Now they are defined with respect to the lowest stateionary state.
	//
	binning = 0; //binning counter
	double* nn;
	double test;
	if (generate_energy_landscape == 1)
	{
		for (int ensemble = 0; ensemble < ensemble_number; ensemble++) // loop over molecules
		{
			GenerateEnergyLandscape(coeff_matrixe, coeff_matrixg, eigenve, eigenvg, TE_matrix_excited, TE_matrix_ground, // transmission probabilities
				attempt_freq_matrix_e, attempt_freq_matrix_g, parabola, energy_difference,
				ensemble, starting_well, binning, v_barrier_height_excited, v_barrier_height_ground, // the latter is used only to report distributions of barrier heights //
				excited_parabola_for_rate, ground_parabola_for_rate,
				r, produced_ensemble, translation, cutoff_energy_difference);
			//
			// one call to GenerateEnergyLandscape produces ONE energy landscape for ONE pigment-protein system, including both ground and excited states. 
			// all low-resolution landscapes are stored in a huge landscape_file matrix.
			// High-resolution landscapes are not saved, only the T-independent tunneling probabilities,attempt frequencies,
			// and coefficients rho are saved. They are used in create_rate_matrix().
			//
			nn = (double*)calloc(point_numbers, sizeof(double));
			for (int i = 0; i < point_numbers; i++)
			{
				test = gsl_matrix_get(ground_parabola_for_rate, ensemble, i);// / 1.98e-23;
				nn[i] = test; // reporting in MHz
			}
			free (nn);

			nn = (double*)calloc(ensemble_number * well_numbers * Hamiltonian_well_numbers, sizeof(double));
			for (int i = 0; i < ensemble_number * well_numbers * Hamiltonian_well_numbers; i++)
			{
				test = gsl_vector_get(coeff_matrixe, i); // projections, for some reason only for lowest group of states
				nn[i] = test;
			}
			
			for (int i = 0; i < ensemble_number * well_numbers * Hamiltonian_well_numbers; i++)
			{
				test = gsl_vector_get(TE_matrix_excited, i); // tunneling probabilities, 10^-25
				nn[i] = test;
			}
			for (int i = 0; i < ensemble_number * well_numbers * Hamiltonian_well_numbers; i++)
			{
				test = gsl_vector_get(attempt_freq_matrix_e, i); // attempt frequencies, slightly less than 10^12
				nn[i] = test;
			}
			free(nn);
		}	// end of for for ensembles
		
	} //end if
	// energy landscapes and T-independent contributions to rates have been generated
	//
	// now that we created the landscapes, we can populate them...
	//
	//How much data will we be storing?
	if (starting_from_equilibrium == 1 || starting_from_uniform_distribution == 1)
		number_of_row_cooling = 1; // in this case the prob_file_cooling matrix will have two rows: one for RT and one for BT
	else
	{ 
	number_of_row_cooling = (int)((initial_temperature-final_temperature)/temperature_step); // number of steps in cooling for which we save probabilities
	}
	// the number of row in the prob_file.txt is related to the integer answer of dividing cooling_time by cooling_time_step 
	//(this is the step numbers which takes the sample to be cooled down)
	//
	gsl_matrix* prob_file_cooling = gsl_matrix_calloc((number_of_row_cooling + 2), number_of_column);
	// the array's size should be number of rows plus 1 because it starts from time zero / initial situation. 
	// it shoud be mentioned that gsl_matrix_alloc (size_t n1, size_t n2) 
	//creates a matrix of size n1 rows by n2 columns

	//next two matrixes are used to keep boltzmann (equilibrium) distributions for room T and low T, 

	gsl_matrix* boltzman_distribution_BT = gsl_matrix_calloc(well_numbers, ensemble_number);
	// this matrix is similar to last_rates or output.txt however its elements (probability in each well) are calculated 
	//according to the Boltzmann distribution at Burn Temperature.
	gsl_matrix* boltzman_distribution_RT = gsl_matrix_calloc(well_numbers, ensemble_number);
	
	//
	gsl_matrix* total_frequency = gsl_matrix_calloc(well_numbers * ensemble_number, 2);
	// this array contains frequencies and respective probabilities, gets created from energy_difference and last row of the respective rate file
	// for the purpose of writing spectra
	// 
	//TO DO? modified transition frequencies based on the lowest eigenstates?
	//

	//these are current rate matrices for current system
	
	gsl_matrix* rate_matrix = gsl_matrix_calloc(well_numbers, well_numbers);
	gsl_matrix* temp_matrix = gsl_matrix_calloc(well_numbers, well_numbers); //backup of the rate matrix
	// after solving the rate equation (SolveRateMatrix), rate_matrix will be changed so Mehdi defined temp_matrix as backup
	
	// The A-s are the probabilities to find the system in a particular well, P-s in our 2015 paper
	A0_ground = gsl_vector_calloc(well_numbers); 
	A0_excited = gsl_vector_calloc(well_numbers);
	
	pre_burn = gsl_matrix_calloc(ensemble_number, well_numbers);// backup of all A0_ground before burn. // matrix as of Nov 26, 2021
	post_burn = gsl_matrix_calloc(ensemble_number, well_numbers);// backup of all A0_ground after burn. 
	//Needed if we try different recovery scenarios with resetting to post-burn situation
	
	int temperature_counter = 0;// this counter is meant to count the number of changes of temperature upon cooling
	
	well_numbers_counter2 = 1; // first (zero-th) item in each row is time
	sum_resonant_probabilities_preburn = 0;

	nn = (double*)calloc(well_numbers, sizeof(double));
	//cooling... BBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBBbb
	for (int ensemble = 0; ensemble < ensemble_number; ensemble++) //for each ensemble
	{
		CalculateBoltzmanDist(ground_parabola_for_rate, boltzman_distribution_RT, ensemble, room_temperature);
		// it calculates the equilibrium ground state distribution at room temperature 

		// *****************************initial condition for the ground state.*****************************
		// It doesn't really matter what the initial distribution is, 
		// since after about 200 microseconds at room temperature the distribution in the ground state 
		//will evolve into Boltzman distribution at room temperature (e.g. ~300K).
		gsl_vector_set_zero(A0_ground);

		gsl_matrix_get_col(A0_ground, boltzman_distribution_RT, ensemble); //boltzmann dist matrixes are well number, ensemble number.
		// ground state probability distribution is set to Boltzmann (equilibrium) distribution for this particular system
		for (int i = 0; i < well_numbers; ++i) // populate zero'th row with initial RT probabilities
		{
			gsl_matrix_set(prob_file_cooling, 0, well_numbers_counter2, gsl_vector_get(A0_ground, i));
			++well_numbers_counter2; // this does not get reset for new molecule
		}

		//*********************************************************************************************************************************************************************

		
		//
		//COOLING******************** TEMPERATURE LOOP *************************
		//
		
		temperature_counter = 0; // zero-th row of prob_file_cooling is populated with room-T probabilities, always
		if (starting_from_equilibrium == 0 && starting_from_uniform_distribution == 0)
		{
			// Then the pre burn distribution is the result of a cooling cycle. The resulting distribution is far from equilibrium.
			gsl_matrix_get_col(A0_ground, boltzman_distribution_RT, ensemble);
			// it sets A0_ground according to boltzmann distribution at room temperature instead of uniform distribution.
			//
			for (double temperature = initial_temperature; temperature >= final_temperature; temperature -= temperature_step)
				// for each temperature, give the system's ground-state probability distribution some time to evolve...
			{
				gsl_matrix_set_zero(rate_matrix); // matrixes need to be reset every time they change,
				gsl_matrix_set_zero(temp_matrix); // or shit gets added up after SolveRateMatrix
				CreateRateMatrix(coeff_matrixg, eigenvg, attempt_freq_matrix_g, TE_matrix_ground, ground_parabola_for_rate, 
					rate_matrix, temp_matrix, temperature, ensemble, v_barrier_height_ground);
				// tunneling and barrier-hopping rates are temperature-dependent, so they have to be recalculated for every temperature
				// cooling is happening in the ground state

				for (int i = 0; i < well_numbers; i++)
				{
					test = gsl_vector_get(A0_ground, i); // probabilities
					nn[i] = test;
				}
				SolveRateMatrix(rate_matrix, ensemble, prob_file_cooling, cooling_time, A0_ground, temperature_counter, temperature);
				// A0_ground are changed from this point on. 
				// need to check what is written in rate_file_cooling
				for (int i = 0; i < well_numbers; i++)
				{
					test = gsl_vector_get(A0_ground, i); // probabilities
					nn[i] = test;
				}
				temperature_counter = temperature_counter + 1;
			}
		}
		else 
		{ // not a rigorous calculation, starting from thermal equilibrium at (low) burn time...
			if (starting_from_equilibrium == 1)
			{
				CalculateBoltzmanDist(ground_parabola_for_rate, boltzman_distribution_BT, ensemble, burn_temperature);
				// need to do it only if we decide to start burning from equilibrium at low T
				// it calculates the equilibrium ground state distribution at burn temperature 
				gsl_matrix_get_col(A0_ground, boltzman_distribution_BT, ensemble);
			}
	
			if (starting_from_uniform_distribution == 1) // not particularly realistic. flash-freeze would freeze boltzmann at RT.
			{
				for (int i = 0; i < well_numbers; i++)
				{
					gsl_vector_set(A0_ground, i, 1 / well_numbers);
				}
			}
		}
		
		nn = (double*)calloc(well_numbers , sizeof(double));
		for (int i = 0; i < well_numbers; i++)
		{
			test = gsl_vector_get(A0_ground, i);// / 1.98e-23;
			nn[i] = test; // reporting probabilites
		}
		//END OF COOLING SAMPLE BEFORE THE OPTICAL EXPERIMENT************************************************************************************
		// save into pre_burn, global variable available to other subroutines:
		gsl_matrix_set_row(pre_burn, ensemble, A0_ground);// A0_ground is now the result of the calculation simulating cooling the sample down

		// preburn matrix has ensemble_number rows and well_numbers columns
		// we are backing it up into pre_burn so we can redo the burn many times with different parameters/conditions

		sum_resonant_probabilities_preburn = sum_resonant_probabilities_preburn + gsl_vector_get(A0_ground, gsl_vector_get(starting_well,ensemble));
		//
		
	} // end of cycle over the molecules.
	free(nn);
	// the pre-burn situation for all molecule+landscape systems is stored in pre-burn matrix containing all A0_ground
	
	gsl_matrix* total_frequency_preburn = gsl_matrix_calloc(well_numbers * ensemble_number, 2);
	gsl_matrix* last_probs_preburn = gsl_matrix_calloc(well_numbers, ensemble_number);// this is not necessary. Mehdi had it for debugging 

	//total_bins = (int)sqrt(well_numbers * ensemble_numbery);
	bin = gsl_vector_calloc(total_bins + 1);
	// as the maximum frequency should be considered in binning, there is a necessity to go one step beyond the range of frequency.
	
	ConvertWellToFrequency(energy_difference, prob_file_cooling, last_probs_preburn, total_frequency_preburn, number_of_row_cooling);
	// this thing takes energy_difference and prob_file type data and creates last_probs type matrix (again "rates" mean probabilities) and 
	// total frequency type matrix, with two rows, one containing all frequencies and another - all probabilities

	BinFrequency(total_frequency_preburn, bin, bin_step, smallest_frequency);
	//this thing uses as input the total_frequency file that contains frequencies and probabilities.
		// it calculates minimal frequency, maximal frequency and the bin step
		// AND actually bins the probabilites and puts the sums of probabilities into the vector bin. 
		// the frequency information apparently gets recreated elsewhere from smallest_frequency and bin step
	//
	
	WriteLastRowOfRateEquation(last_probs_preburn, "output-preburn.txt"); 
	// looks like this thing is saving some probabilities, and not spectra. most likely the final probabilities after cooling

	WriteBinnedFrequency(bin, bin_step, smallest_frequency, "preburnspectrum.txt");// frequency_binned.txt in the previous versions.

	//
	// //
	// Ultra-long evolution of absorption spectrum. In principle, could be modeled using recovery 
	//without burn to do the same thing
	// only makes sense to do if realistic cooling was employed first
	if (starting_from_equilibrium == 0 && starting_from_uniform_distribution == 0)
	{
		for (int ensemble = 0; ensemble < ensemble_number; ensemble++)
		{
			gsl_matrix_set_zero(rate_matrix); // matrixes need to be reset every time they change,
			gsl_matrix_set_zero(temp_matrix); // or shit gets added up after SolveRateMatrix
			CreateRateMatrix(coeff_matrixg, eigenvg, attempt_freq_matrix_g, TE_matrix_ground, ground_parabola_for_rate,
				rate_matrix, temp_matrix, burn_temperature, ensemble, v_barrier_height_ground);
			SolveRateMatrix(rate_matrix, ensemble, prob_file_cooling, whole_long_time, A0_ground, temperature_counter, burn_temperature);
			// whole_long_time = 13 years
		}
		ConvertWellToFrequency(energy_difference, prob_file_cooling, last_probs_preburn, total_frequency_preburn, number_of_row_cooling);
		// this thing takes energy_difference and prob_file type data and creates last_probs type matrix (again "rates" mean probabilities) and 
		// total frequency type matrix, with two rows, one containing all frequencies and another - all probabilities

		BinFrequency(total_frequency_preburn, bin, bin_step, smallest_frequency);
		//this thing uses as input the total_frequency file that contains frequencies and probabilities.
			// it calculates minimal frequency, maximal frequency and the bin step
			// AND actually bins the probabilites and puts the sums of probabilities into the vector bin. 
			// the frequency information apparently gets recreated elsewhere from smallest_frequency and bin step

		WriteBinnedFrequency(bin, bin_step, smallest_frequency, "longtermspectrum.txt");// frequency_binned absorption spectrum for ulta-long time
	}

	// TO DO: free memory!
	gsl_rng_free(r);
	gsl_matrix_free(rate_matrix);
	gsl_matrix_free(temp_matrix);
	gsl_matrix_free(last_probs_preburn);
	gsl_matrix_free(total_frequency_preburn);
	// various items necessary to recreate the rate matrix at different temperatures are not erased.
	 
}


void CMitacsGUIDlg::OnBnClickedCancel()
{
	// TODO: Add your control notification handler code here
	CDialogEx::OnCancel();
}

// The next several ones are checkboxes, that are handled as buttons
void CMitacsGUIDlg::OnBnClickedEquilibrium()
{
	if (Equilibrium.GetCheck() == BST_CHECKED)
	{
		starting_from_equilibrium = 1; //getting here if button was unchecked
		starting_from_uniform_distribution = 0;
			}
	else
	{
		starting_from_equilibrium = 0;
		starting_from_uniform_distribution = 1; //getting here if button was checked
		Uniform.SetCheck(BST_CHECKED);
		UpdateData();
	}
}

void CMitacsGUIDlg::OnBnClickedUniform()
{
	if (Uniform.GetCheck() == BST_CHECKED)
	{
		starting_from_uniform_distribution = 1;//getting here if unchecked before click
		starting_from_equilibrium = 0;
	}
	else
	{
		starting_from_uniform_distribution = 0;
		starting_from_equilibrium = 1;
		Equilibrium.SetCheck(BST_CHECKED);
		UpdateData();
	}
}

void CMitacsGUIDlg::OnBnClickedgenerate()
{
	if (generate.GetCheck() == BST_CHECKED)
	{
		generate_energy_landscape = 1;//getting here if unchecked before click
			}
	else
	{
		generate_energy_landscape = 0;
			}
	}

void CMitacsGUIDlg::OnBnClickedbalance()
{
	if (balance.GetCheck() == BST_CHECKED)
	{
		balancing = 1;//getting here if unchecked before click
	}
	else
	{
		balancing = 0;
	}
}

void CMitacsGUIDlg::OnBnClickedbottom()
{
	if (bottom.GetCheck() == BST_CHECKED)
	{
		bottom_decoupling = 1;//getting here if unchecked before click
	}
	else
	{
		bottom_decoupling = 0;
	}
}

void CMitacsGUIDlg::OnBnClickedbarrier()
{
	if (barrier.GetCheck() == BST_CHECKED)
	{
		barrier_decoupling = 1;//getting here if unchecked before click
	}
	else
	{
		barrier_decoupling = 0;
	}
}

void CMitacsGUIDlg::OnBnClickedanalytical()
{
	if (analytical.GetCheck() == BST_CHECKED)
	{
		analytical_burn_calculation = 1;//getting here if unchecked before click
		recovery_between_two_acts_of_burn = 0;
		recovery_burn.SetCheck(BST_UNCHECKED);
		burn_in_one_step = 0;
		one_step.SetCheck(BST_UNCHECKED);
	}
	else
	{
		analytical_burn_calculation = 0;
	}
}

void CMitacsGUIDlg::OnBnClickedonestep()
{
	if (one_step.GetCheck() == BST_CHECKED)
	{
		burn_in_one_step = 1;//getting here if unchecked before click
		recovery_between_two_acts_of_burn = 0;
		recovery_burn.SetCheck(BST_UNCHECKED);
		analytical_burn_calculation = 0;
		analytical.SetCheck(BST_UNCHECKED);
	}
	else
	{
		burn_in_one_step = 0;
	}
}

void CMitacsGUIDlg::OnBnClickedrecoveryburn()
{
	if (recovery_burn.GetCheck() == BST_CHECKED)
	{
		burn_in_one_step = 0;//getting here if unchecked before click
		recovery_between_two_acts_of_burn = 1;
		analytical_burn_calculation = 0;
		analytical.SetCheck(BST_UNCHECKED);
		one_step.SetCheck(BST_UNCHECKED);
	}
	else
	{
		recovery_between_two_acts_of_burn = 0;
	}
}

void CMitacsGUIDlg::OnBnClickedscan()
{
	if (scan.GetCheck() == BST_CHECKED)
	{
		recovery_correction_during_frequency_scan = 1;//getting here if unchecked before click
	}
	else
	{
		recovery_correction_during_frequency_scan = 0;
	}
}

void CMitacsGUIDlg::OnBnClickeddetails()
{
	if (details.GetCheck() == BST_CHECKED)
	{
		writing_details_of_burning_in_reson_well = 1;//getting here if unchecked before click
	}
	else
	{
		writing_details_of_burning_in_reson_well = 0;
	}
}

// Checkbox region has ended.
//Below we have code for actual buttons.
//
//
//burning a hole
void CMitacsGUIDlg::OnBnClickedBurn() 
// supposed result of running this subroutine is HGK plus final absorption or hole spectrum
// OR some array or matrix from which these can be generated.
//HGK: average (over molecules) population of the resonant well, as a function of time or step number
//hole: populations in all wells for all molecules after the burn
{
	double prob_in_reson_well_after_each_burn_counter, avg_prob_reson_well_counter, sum_over_probability_in_resonant_well;
	double simulation_time;
	double burning_yield_in_resonant_well, hole_burning_yield_summation, test;
	int st_well;
	double* mm;
	sum_over_probability_in_resonant_well = 0;
	hole_burning_yield_summation = 0;
	double hole_burning_yield = 0;
	int HGK_counter = 0; // looks like the index of the step in HGK
	gsl_matrix* rate_matrix_g = gsl_matrix_calloc(well_numbers, well_numbers);
	gsl_matrix* temp_matrix_g = gsl_matrix_calloc(well_numbers, well_numbers); //backup of the rate matrix - ground state
	gsl_matrix* rate_matrix_e = gsl_matrix_calloc(well_numbers, well_numbers);
	gsl_matrix* temp_matrix_e = gsl_matrix_calloc(well_numbers, well_numbers); //backup of the rate matrix
	gsl_vector* A0_excited_backup = gsl_vector_calloc(well_numbers); 
	// backup of the A0_excited vector after evolving from 000010000. Allows to not recalculate excited state evolution
		
	gsl_vector* hole_burning_yield_starting_well = gsl_vector_calloc(ensemble_number);
	// It keeps the probability of leaving starting well after burning / after one step of evolution in the excited state 
	//(which is something smaller than one).
	gsl_vector* pre_burn_at_starting_well = gsl_vector_calloc(ensemble_number);
	// It keeps the probability for starting well in the ground state before burning

	
	//
	gsl_matrix* boltzman_distribution_BT = gsl_matrix_calloc(well_numbers, ensemble_number);
	gsl_matrix* hole_wells = gsl_matrix_calloc(well_numbers, ensemble_number);
	// here we save post_burn and pre_burn differences, with A0 being columns

	gsl_matrix* hole_growth_kinetics;
		
	number_of_HGK_measurement = (int)(number_of_burn / HGK_measurement_frequency);
	//The number of points in the reported HGK curve

	if (burn_in_one_step == 1) // prob_file_burning will contain two rows, pre-burn and post-burn
	{
		hole_growth_kinetics = gsl_matrix_calloc(1, ensemble_number);
		number_of_row_burning = 1;
	}
	// in case that the whole burning is calculated in one step (we are only interested in the final hole), this matrix has only one row
	else
	{
		hole_growth_kinetics = gsl_matrix_calloc(number_of_HGK_measurement+1, ensemble_number);
		number_of_row_burning = number_of_HGK_measurement;
	}
	// matrix storing results of hole growth kinetics modeling, separately for every molecule
	//
	gsl_matrix* prob_file_burning = gsl_matrix_calloc((number_of_row_burning + 1), number_of_column);
	gsl_matrix* prob_file_burning_e = gsl_matrix_calloc(1, number_of_column); // the place for SolvRateMatrix to dump a0_excited...
	// to avoid dumping it into prob_file_burning

	int well_numbers_counter2 = 1;
	sum_resonant_probabilities_postburn = 0;
	prob_in_reson_well_after_each_burn_counter = 0;
	avg_prob_reson_well_counter = 0;

	mm = (double*)calloc(point_numbers, sizeof(double));
	for (int ensemble = 0; ensemble < ensemble_number; ensemble++) // for each molecule
	{
		HGK_counter = 0;

		//starting well
		st_well = gsl_vector_get(starting_well, ensemble); 
		// the most resonant well, can change from molecule to molecule. In fixed-frequency burning it does not change until molecule change
		//
		// Initial conditions for burning in terms of pre-burn distributions A0_ground
		if (starting_from_equilibrium == 1) // recalculates the pre-burn probability distribution, it is no longer the result of realistic cooling.
		{
			CalculateBoltzmanDist(ground_parabola_for_rate, boltzman_distribution_BT, ensemble, burn_temperature);
			// need to do it only if we decide to start burning from equilibrium at low T
			// it calculates the equilibrium ground state distribution at burn temperature 
			gsl_matrix_get_col(A0_ground, boltzman_distribution_BT, ensemble);
		}

		if (starting_from_uniform_distribution == 1) // starting burn from uniform probability distribution. 
			//One will have this situation in case of hyperquenching / flash freezing?
			// or, rather, one should be using room-T Bolzmann distribution here.
		{
			for (int i = 0; i < well_numbers; i++)
			{
				gsl_vector_set(A0_ground, i, 1 / well_numbers);
			}
		}
		//
		
		//if not doing one of special cases, just get A0 from backup (pre-burn / result of cooling process)...
		if (starting_from_equilibrium == 0 && starting_from_uniform_distribution == 0)
		{
			gsl_matrix_get_row(A0_ground, pre_burn, ensemble); // get it from pre-burn
		}
		//
		//
		//
		for (int i = 0; i < well_numbers; i++) // populate zero'th row of prob_file_burning with pre-burn probabilities
		{
			gsl_matrix_set(prob_file_burning, 0, well_numbers_counter2, gsl_vector_get(A0_ground, i));
			++well_numbers_counter2; // this does not get reset for new molecule
		}
		// 
		//Excited state
		gsl_matrix_set_zero(rate_matrix_e); // matrixes need to be reset every time they change,
		gsl_matrix_set_zero(temp_matrix_e); // or shit gets added up after SolveRateMatrix
		CreateRateMatrix(coeff_matrixe, eigenve, attempt_freq_matrix_e, TE_matrix_excited, excited_parabola_for_rate,
			rate_matrix_e, temp_matrix_e, burn_temperature, ensemble, v_barrier_height_excited);
		// tunneling and barrier-hopping rates are temperature-dependent, but burning occurs at fixed temperature, so we only need to do it once.
		////
		for (int i = 0; i < well_numbers; i++)
		{
			test = gsl_matrix_get(rate_matrix_e, 2, i); // i-th row of rate matrix
			mm[i] = test;
		}
		
		if (recovery_between_two_acts_of_burn == 1)
		{    //recovery between two consecutive acts of the photon absorption - in the ground state
			gsl_matrix_set_zero(rate_matrix_g); // matrixes need to be reset every time they change,
			gsl_matrix_set_zero(temp_matrix_g); // or shit gets added up after SolveRateMatrix
			CreateRateMatrix(coeff_matrixg, eigenvg, attempt_freq_matrix_g, TE_matrix_ground, ground_parabola_for_rate,
				rate_matrix_g, temp_matrix_g, burn_temperature, ensemble, v_barrier_height_ground);
			// recovery while burning is happening at fixed temperature, same as burn temperature, so need to do it only once
			for (int i = 0; i < well_numbers; i++)
			{
				test = gsl_matrix_get(rate_matrix_g, 2, i); // i-th row of rate matrix
				mm[i] = test;
			}
		}
		//
		//
		
		if (burn_in_one_step == 0) // doing rigorous step by step burning, possibly including the recovery while burning
		{
			// modification for faster calculations, December 2022: every evolution in the exited state starts in exactly same condition. 
			// and takes the same time, so one does not need to re-solve exactly the same excited state matrix all over again.
			
			for (int burn_counter = 0; burn_counter < number_of_burn; burn_counter++) // do as many times as there are individual acts of burn
			{
				if (burn_counter == 0) // do it only the first time, and then keep evolved A0_excited_backup
					// for every act of burn this step is exactly the same
				{
					gsl_vector_set_zero(A0_excited); // all probabilities in the excited state are set to zero
					gsl_vector_set(A0_excited, st_well, 1); // starting probability in the resonant well of the excited state is set to one
					gsl_matrix_memcpy(rate_matrix_e, temp_matrix_e); // destination, source
					SolveRateMatrix(rate_matrix_e, ensemble, prob_file_burning_e, burn_time, A0_excited, burn_counter+1, -1000);
					// -1000 is pseudo-temperature, used by Mehdi to see at which instance of solveRateMatrix it crashed
					// burn_time is excited state lifetime
					// probabilities are allowed to evolve in the excited state; NB: A0_excited are changed from this point on. 
					// NB: prob_file_burning_e will contain A0-excited and we do not care, wedo not use it, just dump stuff there
					//
					burning_yield_in_resonant_well = 1 - gsl_vector_get(A0_excited, st_well);
					hole_burning_yield_summation += burning_yield_in_resonant_well;
					// the difference between the probability to be in resonant well before and after burn can be interpreted as yield.
					gsl_vector_memcpy(A0_excited_backup, A0_excited);
				}
				else gsl_vector_memcpy(A0_excited, A0_excited_backup); // if it is not the first step of burn
				
				//and the following is done on every step
				gsl_vector_scale(A0_excited, gsl_vector_get(A0_ground, st_well)); // essentally multiplying all parts of the modified excited 
				//state probability vector by the probability to be in the resonant well of the ground state

				gsl_vector_set(A0_ground, st_well, 0); // erase preburn probability in the resonant well of the ground state
				gsl_vector_add(A0_ground, A0_excited); // update the probabilities in the ground state

				// TO DO: check if new probabilities add up to one
				//
				if (recovery_between_two_acts_of_burn == 1)
				{    //recovery between two consecutive acts of the photon absorption.
					simulation_time = time_interval_between_two_consecutive_burns * 1e6; 
					// in microseconds... apparently this is what solve rate matrix wants
					// this time is calculated in seconds earlier, however time scale in this program elsewhere is in microsecond.
					//
					gsl_matrix_memcpy(rate_matrix_g, temp_matrix_g); // reset rate matrix to pre-evolution state
					SolveRateMatrix(rate_matrix_g, ensemble, prob_file_burning, simulation_time, A0_ground, burn_counter+1, -110);
					// A0_ground are changed from this point on. 
					// the rate_matrix and A0_ground are changed. NOTE: we want to calculate the distribution in each time step before 
					//the next photon is absorbed by the pigment. By -110 as a temperature, if something happens by mistake in this function 
					// Mehdi could catch it.
					// if this is done, the prof_file_burning must contain proper numbers corresponding to A0_ground
				}
				else // populate respective row with contents of A0_ground after burn
				{
					for (int i = 0; i < well_numbers; i++) 
					{
						gsl_matrix_set(prob_file_burning, burn_counter + 1, ensemble*well_numbers+i+1, gsl_vector_get(A0_ground, i));
						
					}
				}
				//
				//sum_over_probability_in_resonant_well += gsl_vector_get(A0_excited, st_well); // sum over molecules?
				//
				// HOLE GROWTH KINETICS
				//
				if ((burn_counter) % HGK_measurement_frequency == 0 && number_of_HGK_measurement != 1)
				{		
					//save results into H_G_K array every now and then, not after every act of burn...
					gsl_matrix_set(hole_growth_kinetics, HGK_counter, ensemble, gsl_vector_get(A0_ground, st_well));
					++HGK_counter;
				}
				if (number_of_HGK_measurement == 1)
					// in case that we want to measure the hole once only, the last result should be saved for burning not the first one.
					gsl_matrix_set(hole_growth_kinetics, 0, ensemble, gsl_vector_get(A0_ground, st_well));
				// so the row of the h_g_k matrix is HGK for one molecule
				//
			}// the end of burn steps loop.
			 
			//
			//gsl_matrix_set(averaged_probability_in_resonant_well, ensemble, avg_prob_reson_well_counter, sum_over_probability_in_resonant_well / number_of_burn);
			//Not sure what this is... need to debug and see how variables change
			
			//gsl_vector_set(hole_burning_yield_starting_well, ensemble, 1 - gsl_vector_get(A0_excited, st_well));
			gsl_vector_set(pre_burn_at_starting_well, ensemble, gsl_matrix_get(pre_burn, ensemble, st_well));
			//
			// anyway, the end results of any kind of burn, regardless of if we care for HGK, should be the set  
			gsl_matrix_set_row(post_burn, ensemble, A0_ground);
			gsl_matrix_memcpy(post_burn_for_recovery, post_burn); // post_burn will change after SetHoleProbability
		} // end if rigorous burn modeling
		//
		//
	// we are still inside the cycle over molecules
		//
		//
		else // "analytical burn", no recovery while burning, but still multi-step, so HGK is possible
		{
			// calculating the burn process via analytical equations (see supplemental of the paper) 
			//In case that recovery between two consecutive acts of burn is on, then analytical model is not valid.
			if (number_of_burn == 0)
			{
				//cout << "The system should experience at least one act of NPHB" << endl;
				throw;
			}
			double probability_excited_after_resonant, probability_ground_before_resonant, probability_excited_after_nonresonant, probability_ground_before_nonresonant;
			double probability, z, series_answer;
			int counter_burn_number, summation_burn_number;
			gsl_matrix* prob_in_reson_well_after_each_act_burn;
			prob_in_reson_well_after_each_act_burn = gsl_matrix_calloc(ensemble_number, well_numbers * (number_of_burn + 1) );
			gsl_matrix* averaged_probability_in_resonant_well;
			gsl_matrix* probabilities_in_each_frequency_step;
			//
			//starts the same way as the algorithm above.
			//
			gsl_vector_set_zero(A0_excited);
			gsl_vector_set(A0_excited, st_well, 1); // excited state probabilities start off with 1 in the resonant well and 
			//zero everywhere else
			//gsl_vector_set (A0_excited, starting_well, gsl_vector_get (A0_ground, starting_well));
			//
			gsl_matrix_memcpy(rate_matrix_e, temp_matrix_e); // destination, source
			SolveRateMatrix(rate_matrix_e, ensemble, prob_file_burning_e, burn_time, A0_excited, number_of_row_burning, -1000);
			// the rate_matrix and A0_excited are changed from this point on. 
			//
			// after evolving for the [excited state lifetime] the AO_excited got modified and probability in the resonant well is lower than one
			//
			burning_yield_in_resonant_well = 1 - gsl_vector_get(A0_excited, st_well);
			// the difference between the probability at resonant well before (1) and after burn can be interpreted as effective hole burning yield.
			//
			hole_burning_yield_summation += burning_yield_in_resonant_well; // sum for all molecules?
			// most likely we are summing it up to later calculate the averge yield for resonant well over all systems ...
			//
			gsl_vector_set(hole_burning_yield_starting_well, ensemble, burning_yield_in_resonant_well);
			// saving  HB yields for every molecule
			//
			probability_excited_after_resonant = gsl_vector_get(A0_excited, st_well);
			probability_ground_before_resonant = gsl_vector_get(A0_ground, st_well);
			//
			gsl_vector* A0_excited_temp = gsl_vector_calloc(well_numbers);//it keeps the A0_excited temporarily to use in the following loop.
			gsl_vector_memcpy(A0_excited_temp, A0_excited);
			//
			// in analytical method, probability in wells after only the first burn is considered, and these numbers are reused...
			//
			//
			for (int burn_number = 1; burn_number <= number_of_burn; burn_number++)
				//burn_number is number of the current step of burning
			{
				//++counter_burn_number; // which burn step we are currently on
				//summation_burn_number += burn_number; 
				// this is a weird quantity as it is sum like 1+2+3+4+5... not sure what is the meaning of it
				// need to check what is normalized using this quantity and why...
				//
				for (int k = 0; k < well_numbers; k++)
				{
					//resonant well:
					if (k == st_well)
					{
						probability = pow(probability_excited_after_resonant, burn_number) * probability_ground_before_resonant;
						gsl_vector_set(A0_excited, k, probability);// k is starting_well
						continue;// next k will run.
					}
					// and for non-resonant wells:
					probability_excited_after_nonresonant = gsl_vector_get(A0_excited, k);
					probability_ground_before_nonresonant = gsl_matrix_get(pre_burn, ensemble, k);
					z = probability_excited_after_resonant;
					series_answer = (-1 + pow(z, burn_number)) / (-1 + z);
					probability = probability_excited_after_nonresonant * probability_ground_before_resonant * series_answer + probability_ground_before_nonresonant;
					gsl_vector_set(A0_excited, k, probability);
				} // end of for well numbers
				//
				//
				sum_over_probability_in_resonant_well += gsl_vector_get(A0_excited, st_well);
				// this is the sum of probabilities to remain in the resonant well...
				//
				if (writing_details_of_burning_in_reson_well == 1 && number_of_burn < 100)
				{
					gsl_matrix_set(prob_in_reson_well_after_each_act_burn, ensemble, prob_in_reson_well_after_each_burn_counter, gsl_vector_get(A0_excited, st_well));
					++prob_in_reson_well_after_each_burn_counter;
				}
				gsl_vector_memcpy(A0_ground, A0_excited);
				if ((burn_number - 1) % HGK_measurement_frequency == 0 && number_of_HGK_measurement != 1)
				{
					gsl_matrix_set(hole_growth_kinetics, HGK_counter, ensemble, gsl_vector_get(A0_ground, st_well));
					// burn_number is actually the burn_counter in the versions 20-23. I changed its name since version 19, 
					//I used burn_number to calculate the post-burn analytically. It should be mention that the well_numbers of hole_growth_kinetics is from 0 to number_of_burn-1. Since the loop should start from 1 to be consistent with version 19, burn_number -1 should be used for the index.
					++HGK_counter;
				}
				if (number_of_HGK_measurement == 1) gsl_matrix_set(hole_growth_kinetics, 0, ensemble, gsl_vector_get(A0_ground, st_well));
				// in case that we want to measure the HGK once only, the last result should be saved for burning not the first one.

				// burn_number is actually the burn_counter in the versions 20-23. I changed its name since version 19, 
				//I used burn_number to calculate the post-burn analytically. It should be mention that the well_numbers of hole_growth_kinetics is from 0 to number_of_burn-1. Since the loop should start from 1 to be consistent with version 19, burn_number -1 should be used for the index.
				gsl_vector_memcpy(A0_excited, A0_excited_temp);
			} // end of for many buring steps
			//
			//
			gsl_vector_free(A0_excited_temp);
			gsl_vector_set(pre_burn_at_starting_well, ensemble, gsl_matrix_get(pre_burn, ensemble, st_well));
			gsl_matrix_set_row(post_burn, ensemble, A0_ground);
						
			++avg_prob_reson_well_counter;
		} // else if analytical burn

		sum_resonant_probabilities_postburn = sum_resonant_probabilities_postburn + gsl_vector_get(pre_burn_at_starting_well, ensemble);
	} //end of for, molecules loop
	free(mm);
	gsl_matrix_memcpy(post_burn_for_recovery, post_burn); // post_burn will change after calculating holes; so we are backing it up.
	
	//Replacement for SetHoleProbability(pre_burn, post_burn, hole_wells, ensemble):
	gsl_matrix_sub(post_burn, pre_burn);// post_burn = post_burn - pre_burn
	gsl_matrix_transpose_memcpy(hole_wells, post_burn); // has differences of pre- and post-burn A0-s as columns
	// now hole_wells is a matrix too
	//
	gsl_matrix_memcpy(post_burn, post_burn_for_recovery); // just to be safe, recover post-burn
														  
	// saving the results	

	gsl_matrix* total_frequency_postburn = gsl_matrix_calloc(well_numbers * ensemble_number, 2);
	gsl_matrix* last_probs_preburn = gsl_matrix_calloc(well_numbers, ensemble_number);// this is not necessary. Mehdi had it for debugging 

	//total_bins = (int)sqrt(well_numbers * ensemble_numbery); does not change since cooling
	//bin = gsl_vector_calloc(total_bins + 1);
	
	ConvertWellToFrequency(energy_difference, prob_file_burning, last_probs_preburn, total_frequency_postburn, number_of_row_burning);
	// this thing takes energy_difference and prob_file type data and creates last_probs type matrix and 
	// total frequency type matrix, with two rows, one containing all frequencies and another - all final probabilities

	BinFrequency(total_frequency_postburn, bin, bin_step, smallest_frequency);
	//this thing uses as input the total_frequency file that contains frequencies and probabilities.
		// it calculates minimal frequency, maximal frequency and the bin step
		// AND actually bins the probabilites and puts the sums of probabilities into the vector bin. 
		// the frequency information apparently gets recreated elsewhere from smallest_frequency and bin step
	//

	WriteBinnedFrequency(bin, bin_step, smallest_frequency, "postburnspectrum.txt");// frequency_binned.txt in the previous versions.

	WriteLastRowOfRateEquation(last_probs_preburn, "output-postburn.txt");
	// looks like this thing is saving some probabilities, and not spectra. most likely the final probabilities after cooling

	WriteLastRowOfRateEquation(hole_wells, "hole-well probability.txt");

	WriteBurningYieldvsStartingWell(starting_well, hole_burning_yield_starting_well, pre_burn_at_starting_well, energy_difference);
	// the number of initial well is x axis and their probabilities are in the y-axis.


	// hole_growth_kinetics (burn step, ensemble)
	WriteHoleGrowthKinetics(hole_growth_kinetics, pre_burn_at_starting_well);
	//Hole_growth_kinetics is a matrix with number of burn steps rows and ensemble coulumns, and it contains A0 probabilities for resonant well
	// times and doses are calculated based on number and duration of steps in the burning process.

	gsl_matrix_free(total_frequency_postburn);
	gsl_matrix_free(rate_matrix_g);
	gsl_matrix_free(rate_matrix_e);
	gsl_matrix_free(temp_matrix_g);
	gsl_matrix_free(temp_matrix_e);
	gsl_matrix_free(last_probs_preburn);
	gsl_matrix_free(hole_growth_kinetics);
	gsl_matrix_free(hole_wells);
	gsl_vector_free(A0_excited_backup);
	gsl_vector_free(hole_burning_yield_starting_well);
	gsl_vector_free(pre_burn_at_starting_well);

}

// //
// //
// //
//fixed-temperature recovery
void CMitacsGUIDlg::OnBnClickedrecovery()
{
	//whole_recovery_time = whole_recovery_time * 60 * 1000000; //in microseconds
	//number_of_recovery_measurement = _ttoi(recovery_points_text);
	// it means that during the whole_recovery_time, we will stop the evolution of probabilities and calculate the spectum number_of_recovery_measurement (e.g. 100) times.
	//recovery_time = whole_recovery_time / number_of_recovery_measurement; // will be in microseconds
	int st_well;
	gsl_matrix* rate_matrix_g = gsl_matrix_calloc(well_numbers, well_numbers);
	gsl_matrix* temp_matrix_g = gsl_matrix_calloc(well_numbers, well_numbers); //backup of the rate matrix - ground state

	number_of_row_recovery = number_of_recovery_measurement; //recovery_time/recovery_time_step;

	double b, c;
	//
	gsl_matrix* prob_file_recovery = gsl_matrix_calloc((number_of_row_recovery + 1), number_of_column); 

	gsl_matrix* recovery_kinetics = gsl_matrix_calloc(number_of_recovery_measurement, ensemble_number); // A_ground in resonant well for every step and every molecule

	gsl_vector* sum_recovery_kinetics = gsl_vector_calloc(number_of_recovery_measurement); // sum of all A0_ground at given step

	//gsl_matrix_memcpy(post_burn, post_burn_for_recovery); // recovery starts from the end of previous recovery. Can imitate thermocycling step by step...

	
	//
	for (int ensemble = 0; ensemble < ensemble_number; ensemble++) // for each molecule
	{
		gsl_matrix_get_row(A0_ground, post_burn, ensemble); // pre-recovery state of probabilities for given molecule
		st_well = gsl_vector_get(starting_well, ensemble); // the most resonant well, can change from molecule to molecule. 
		//
		gsl_matrix_set_zero(rate_matrix_g); // matrixes need to be reset every time they change,
		gsl_matrix_set_zero(temp_matrix_g); // or shit gets added up after SolveRateMatrix
		CreateRateMatrix(coeff_matrixg, eigenvg, attempt_freq_matrix_g, TE_matrix_ground, ground_parabola_for_rate,
			rate_matrix_g, temp_matrix_g, burn_temperature, ensemble, v_barrier_height_ground);
		// fixed T recovery is happening at fixed temperature, same as burn temperature
		
		for (int recovery_counter = 1; recovery_counter < number_of_recovery_measurement; recovery_counter++) // do as many times as there are individual acts of burn
		{
			gsl_matrix_memcpy(rate_matrix_g, temp_matrix_g); // reset rate matrix to pre-evolution state
			SolveRateMatrix(rate_matrix_g, ensemble, prob_file_recovery, recovery_time, A0_ground, recovery_counter, -110);
			// A0_ground are changed from this point on. 
			// the rate_matrix and A0_ground are changed. NOTE: we want to calculate the distribution in each time step before 
			b = gsl_vector_get(A0_ground, st_well);
			gsl_matrix_set(recovery_kinetics, recovery_counter, ensemble, b);
		}// the end of recovery steps loop.
		//
		// anyway, the end results of any kind of process, should be saved 
		gsl_matrix_set_row(post_burn, ensemble, A0_ground);
	} // end of loop over molecules

	
	for (int recovery_counter = 1; recovery_counter < number_of_recovery_measurement; recovery_counter++)
		// do as many times as there steps of recovery
	{
		c = 0; // sum of all a0_ground for given step
		for (int ensemble = 0; ensemble < ensemble_number; ensemble++) // for each molecule
		{
			b = gsl_matrix_get(recovery_kinetics, recovery_counter, ensemble);
			c = c + b;
		}// 
		gsl_vector_set(sum_recovery_kinetics, recovery_counter, c);
	}// the end of recovery steps loop.
	//
	gsl_matrix* total_frequency_postburn = gsl_matrix_calloc(well_numbers * ensemble_number, 2);
	gsl_matrix* last_probs_recovery = gsl_matrix_calloc(well_numbers, ensemble_number);// this is not necessary. Mehdi had it for debugging 

	//total_bins = (int)sqrt(well_numbers * ensemble_numbery); does not change since cooling
	//bin = gsl_vector_calloc(total_bins + 1);

	ConvertWellToFrequency(energy_difference, prob_file_recovery, last_probs_recovery, total_frequency_postburn, number_of_row_recovery);
	// this thing takes energy_difference and prob_file type data and creates last_probs type matrix and 
	// total frequency type matrix, with two rows, one containing all frequencies and another - all final probabilities

	BinFrequency(total_frequency_postburn, bin, bin_step, smallest_frequency);
	//this thing uses as input the total_frequency file that contains frequencies and probabilities.
		// it calculates minimal frequency, maximal frequency and the bin step
		// AND actually bins the probabilites and puts the sums of probabilities into the vector bin. 
		// the frequency information apparently gets recreated elsewhere from smallest_frequency and bin step
	//

	WriteBinnedFrequency(bin, bin_step, smallest_frequency, "postrecoveryspectrum.txt");

	// do not rewrite postburn_for_recovery. This way we will be able to add more recovery steps.

	//Next part is modified hole growth kinetics writer originally from write_on_files.cpp, line ~900
	fstream file;
	file.open("Recovery-time.txt", ios::out | ios::app);
	if (!file) {
		cout << "Error opening the file Recovery-time.txt";
	}
	file << 0 << "\t" << 1 <<  endl;
	//first point in recovery curve is time=0, hole =1

	double time, norm;
	
	for (int i = 1; i < number_of_HGK_measurement; i++)
		{
		time = recovery_time * i;
		norm = (sum_resonant_probabilities_preburn-gsl_vector_get(sum_recovery_kinetics, i))/(sum_resonant_probabilities_preburn-sum_resonant_probabilities_postburn);
		// normalized recovery
		file << time << "\t" << norm  << endl;
		}
	
	file.close();
	
	gsl_matrix_free(total_frequency_postburn);
	gsl_matrix_free(rate_matrix_g);
	gsl_matrix_free(temp_matrix_g);
	gsl_matrix_free(last_probs_recovery);
	gsl_vector_free(sum_recovery_kinetics);
	
}

//thermocycling
void CMitacsGUIDlg::OnBnClickedThermocycling()
{
	// first, we assume the existence of a file T(t), containing a two-column array.  is minutes, T in Kelvin
	// let's also require that the first time is zero and first temperature is burn temperature

	//the following is copied from the internet - reading from file
	string line;
	string temp = "";
	ifstream inFile("profile.txt");
	int a = 0;
	double b, c;
	double profile[300][2]; // allowing up to 300 entries in the T(t) file
	// open the file stream
	// check if opening a file failed
	if (inFile.fail()) {
		//cerr << "Error opeing a file" << endl;
		inFile.close();
		exit(1);
	}
	while (getline(inFile, line))
	{
		int b = 0; // column index
		for (int i = 0; i < line.size(); i++)
		{ // for each character in rowstring
			if (!isblank(line[i])) { // if it is not blank, do this
				string d(1, line[i]); // convert character to string
				temp.append(d); // append the two strings
			}
			else {
				profile[a][b] = stod(temp);  // convert string to double
				temp = ""; // reset the capture
				b++; // increment b cause we have a new number
			}
		}
		profile[a][b] = stod(temp);
		temp = "";
		a++; // onto next row
	}
	inFile.close(); // close the file stream
	//part copied from the internet ends. TO DO: May want to make a separate subroutine that reads file and returns an array.

	number_of_recovery_measurement = a;
	number_of_row_thermocycling = a;
	//
	int st_well;
	gsl_matrix* rate_matrix_g = gsl_matrix_calloc(well_numbers, well_numbers); // recovery is occuring via ground-statre processes only
	gsl_matrix* temp_matrix_g = gsl_matrix_calloc(well_numbers, well_numbers); //backup of the rate matrix - ground state
	//
	gsl_matrix* prob_file_recovery = gsl_matrix_calloc((number_of_row_thermocycling + 1), number_of_column);

	gsl_matrix* recovery_kinetics = gsl_matrix_calloc(number_of_recovery_measurement, ensemble_number);

	gsl_vector* sum_recovery_kinetics = gsl_vector_calloc(number_of_recovery_measurement); // sum of all A0_ground at given step

	gsl_matrix_memcpy(post_burn, post_burn_for_recovery); // thermocycling always starts from post-burn state. This allows to test multiple thermocycling scenarios after the same burn

	for (int ensemble = 0; ensemble < ensemble_number; ensemble++) // for each molecule
	{
		gsl_matrix_get_row(A0_ground, post_burn, ensemble); // pre-recovery state of probabilities for given molecule
		st_well = gsl_vector_get(starting_well, ensemble); // the most resonant well, can change from molecule to molecule. We assume that we start from here.
		//
		c = 0;
		for (int recovery_counter = 1; recovery_counter < number_of_recovery_measurement; recovery_counter++) // do as many times as there are individual steps in thermocycling
		{
			gsl_matrix_set_zero(rate_matrix_g); // matrixes need to be reset every time they change,
			gsl_matrix_set_zero(temp_matrix_g); // or shit gets added up after SolveRateMatrix
			CreateRateMatrix(coeff_matrixg, eigenvg, attempt_freq_matrix_g, TE_matrix_ground, ground_parabola_for_rate,
				rate_matrix_g, temp_matrix_g, profile[1][recovery_counter], ensemble, v_barrier_height_ground);
			// temperature is changing at every step, so rate matrix is changing at every step
			// //temperature is profile[1][recovery_counter]
			// 
	//TO DO -- calculate time at given temperature / interpolate?
			recovery_time = (profile[0][recovery_counter] - profile[0][recovery_counter - 1]) * 60 * 1000000; // time, the first time in the profile should be zero; converted to microsec
			if (recovery_time < 0) break; // extra safety in case the program tries to process some zero time / zero T line
			SolveRateMatrix(rate_matrix_g, ensemble, prob_file_recovery, recovery_time, A0_ground, number_of_row_thermocycling, -110);
			// A0_ground are changed from this point on. 
			// the rate_matrix and A0_ground are changed. 
			b = gsl_vector_get(A0_ground, st_well);
			gsl_matrix_set(recovery_kinetics, recovery_counter, ensemble, b);
			
		}// the end of recovery steps loop.
		//
		// anyway, the end results of any kind of process, should be recorded 
		gsl_matrix_set_row(post_burn, ensemble, A0_ground); // saving all probability distributions, for the purpose of making spectra

	} // end of loop over molecules
	//
	for (int recovery_counter = 1; recovery_counter < number_of_recovery_measurement; recovery_counter++)
		// do as many times as there steps of recovery
	{
		c = 0; // sum of all a0_ground for given step
		for (int ensemble = 0; ensemble < ensemble_number; ensemble++) // for each molecule
		{
			b = gsl_matrix_get(recovery_kinetics, recovery_counter, ensemble);
			c = c + b;
		}// 
		gsl_vector_set(sum_recovery_kinetics, recovery_counter, c);
	}// the end of recovery steps loop.
	// 
	gsl_matrix* total_frequency_postburn = gsl_matrix_calloc(well_numbers * ensemble_number, 2);
	gsl_matrix* last_probs_recovery = gsl_matrix_calloc(well_numbers, ensemble_number);// this is not necessary. Mehdi had it for debugging 

	//total_bins = (int)sqrt(well_numbers * ensemble_numbery); does not change since cooling
	//bin = gsl_vector_calloc(total_bins + 1);

	ConvertWellToFrequency(energy_difference, prob_file_recovery, last_probs_recovery, total_frequency_postburn, number_of_row_thermocycling);
	// this thing takes energy_difference and prob_file type data and creates last_probs type matrix and 
	// total frequency type matrix, with two rows, one containing all frequencies and another - all final probabilities

	BinFrequency(total_frequency_postburn, bin, bin_step, smallest_frequency);
	//this thing uses as input the total_frequency file that contains frequencies and probabilities.
		// it calculates minimal frequency, maximal frequency and the bin step
		// AND actually bins the probabilites and puts the sums of probabilities into the vector bin. 
		// the frequency information apparently gets recreated elsewhere from smallest_frequency and bin step
	//

	WriteBinnedFrequency(bin, bin_step, smallest_frequency, "postthermocyclingspectrum.txt");

	// do not rewrite postburn_for_recovery. This way we will be able to add more recovery steps.
	// 
	//Next part is modified hole growth kinetics writer originally from write_on_files.cpp, line ~900
	fstream file;
	file.open("Thermocycling-time-T.txt", ios::out | ios::app);
	if (!file) {
		cout << "Error opening the file Thermocycling-time-T.txt";
	}
	file << 0 << "\t" << 1 << "\t" << profile[0][0]<<endl;
	//first point in recovery curve is time=0, hole =1, third line is the first temperature in the profile

	double time, norm;

	for (int i = 1; i < number_of_HGK_measurement; i++)
	{
		time = profile[0][i];
		norm = (sum_resonant_probabilities_preburn - gsl_vector_get(sum_recovery_kinetics, i)) / (sum_resonant_probabilities_preburn - sum_resonant_probabilities_postburn);
		// normalized recovery
		file << time << "\t" << norm <<  "\t" << profile[1][i] << endl;
	}

	file.close();

	gsl_matrix_free(total_frequency_postburn);
	gsl_matrix_free(rate_matrix_g);
	gsl_matrix_free(temp_matrix_g);
	gsl_matrix_free(last_probs_recovery);
	gsl_vector_free(sum_recovery_kinetics);

}



// //
//single molecule scan, modified by Jing in 2023. 
//Unlike Mehdi's quasi-SMS, it has discrete jumps. No evolution based on rate matrix
// instead we calculate HB yields from the rates properly picked from the rate matrix and compare the HB yield with random number to see if jump occurred.
void CMitacsGUIDlg::OnBnClickedSmsscan() 
// This should be run with small number of molecules, 
// as averaging single-molecule level information is not necessary here
		//
{
	double current_frequency, resonant_frequency; // running / constantly changing frequency and current resonant frequency,  
	// that may change as a result of a spectral jump
	double total_time; // total time counted from the very beginning of the first scan, to produce long trajectories.
	int st_well; // can re-use for current resonant well
	double hole_burning_yield_left = 0; // HB yields for light process
	double hole_burning_yield_right = 0;
	double dark_yield_left = 0; // analog of HB yields for dark process
	double dark_yield_right = 0;
	gsl_matrix* rate_matrix_g = gsl_matrix_calloc(well_numbers, well_numbers);
	gsl_matrix* temp_matrix_g = gsl_matrix_calloc(well_numbers, well_numbers);
	gsl_matrix* rate_matrix_e = gsl_matrix_calloc(well_numbers, well_numbers);
	gsl_matrix* temp_matrix_e = gsl_matrix_calloc(well_numbers, well_numbers);
	gsl_matrix* boltzman_distribution_BT = gsl_matrix_calloc(well_numbers, ensemble_number);
	
	
	// how many times in a row we spend in resonance while scanning
	//number_of_burn = (cutoff_energy_difference/ time_interval_between_two_consecutive_burns ) /( (maximum_frequency - minimum_frequency) / scan_time);
	// may not need this if at every step we check if we are in resonance or not.
	int number_of_time_intervals = (int)(scan_time / time_interval_between_two_consecutive_burns);

	CheckFrequencyRange(energy_difference); // modified in 2022, now checks the whole set of transition frequencies

	for (int ensemble = 0; ensemble < ensemble_number; ensemble++) // loop over molecules
	{
	CreateRateMatrix(coeff_matrixg, eigenvg, attempt_freq_matrix_g, TE_matrix_ground, ground_parabola_for_rate,
			rate_matrix_g, temp_matrix_g, burn_temperature, ensemble, v_barrier_height_ground);
		// we assume that scanning occurs at the same temperature, so both burning and recovery occur at the same T
		
		CreateRateMatrix(coeff_matrixe, eigenve, attempt_freq_matrix_e, TE_matrix_excited, excited_parabola_for_rate,
			rate_matrix_e, temp_matrix_e, burn_temperature, ensemble, v_barrier_height_excited);
		// Rate matrix in the excited state - burning
		// NB, for making the program faster: the temperature is NOT changing 
		// //so these rate matrixes could be calculated only once for each molecule
		//
		// //
		// //
		// Various possible starting conditions:
		if (starting_from_equilibrium == 1) // recalculates the pre-burn probability distribution, it is no longer the result of realistic cooling.
		{
			CalculateBoltzmanDist(ground_parabola_for_rate, boltzman_distribution_BT, ensemble, burn_temperature);
			// need to do it only if we decide to start experiemnt from equilibrium at low T
			// it calculates the equilibrium ground state distribution at burn temperature 
			gsl_matrix_get_col(A0_ground, boltzman_distribution_BT, ensemble);
		}
		if (starting_from_uniform_distribution == 1) // starting burn from uniform probability distribution. One will have this situation in case of hyperquenching / flash freezing
			// or, rather, one should be using room-T boltzmann distribution here.
		{
			for (int i = 0; i < well_numbers; i++)
			{
				gsl_vector_set(A0_ground, i, 1 / well_numbers);
			}
		}
		//
		// otherwise uses the results of cooling as a starting point
		//if not doing one of special cases, just get A0 from backup (pre-burn)...
		if (starting_from_equilibrium == 0 && starting_from_uniform_distribution == 0)
		{
			gsl_matrix_get_row(A0_ground, pre_burn, ensemble);
		}
		//
		// 
		// TO DO: Here goes the code where we somehow decide which well is the starting well in the very beginning, before the first scan. 
		// May be completely random, or may be based on one of the above starting population distributions that are now contained in A0_ground.
		// 
		for (int scan_counter = 0; scan_counter < number_of_scan; scan_counter++) // loop over number of scans
		{ 
			// set current frequency to starting frequency of the scan
			// //
			// loop over as many time steps as there are in one scan time (time of scan divided by time between two absorbed photons)
			for (int i = 0; i < number_of_time_intervals; i++)
			{
				// calculate current frequency
				// 
				// check if we are in resonance with the well where the molecule currenlty resides.

				if (recovery_correction_during_frequency_scan == 1 ) // if we take into account dark processes
				{
					//TO DO: code for "tossing a coin" for dark process
					// // use data from the ground-state rate matrix and time_interval_between_two_consecutive_burns to calculate probability of jump in one step (analog of the HB yield)
					//remember that if we have multiple wells, we need to test for left and right jumps
					// if necessary (if jump occurred), reset current resonant well and resonant frequency
					// flag=1 (0 is no jump, 1 is dark jump, 2 is light jump)
				} // endif
				//
				
				// Light processes
				if (recovery_between_two_acts_of_burn == 1) // change this IF to checking for if we are in resonance
				{
					//TO DO: code for "tossing a coin" for light process
					// use data from the excited-state rate matrix and excited state lifetime to calculate probability of jump in one step (HB yield)
					// 
					// remember that if we have multiple wells, we need to test for left and right jumps
					// if necessary, reset resonant well and resonant frequency
					// flag=2 (0 is no jump, 1 is dark jump, 2 is light jump)
				} 
				
				// // save total time counted from the beginning of the first scan, current frequency, resonant well and resonant frequency as well as the "flag" into some matrix or into some file.
				// By "flag" I mean let's say, 0 for no jump, 1 for dark jump and 2 for light-induced jump.
				// Looks like this is all the trajectory information that we need.
				// Can analyze it internally or externally.
				// //

			}// the end of multi-step varying-frequency loop
			//
			
		}	//end of many scan loop
		//
		
	
		
	} // end of loop over molecules (ensemble)
	//TO DO: saving some sort of data / analyzing the trajectories...
}
//

