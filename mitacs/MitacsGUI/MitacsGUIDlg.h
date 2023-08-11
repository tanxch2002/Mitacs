
// MitacsGUIDlg.h: 头文件
//

#pragma once


// CMitacsGUIDlg class
class CMitacsGUIDlg : public CDialogEx
{
// 构造
public:
	CMitacsGUIDlg(CWnd* pParent = nullptr);	// does not have a parent

// 对话框数据
#ifdef AFX_DESIGN_TIME
	enum { IDD = IDD_MITACSGUI_DIALOG };
#endif

	protected:
	virtual void DoDataExchange(CDataExchange* pDX);	// DDX/DDV 支持


// 实现
protected:
	HICON m_hIcon;

	// 生成的消息映射函数
	virtual BOOL OnInitDialog();
	afx_msg void OnSysCommand(UINT nID, LPARAM lParam);
	afx_msg void OnPaint();
	afx_msg HCURSOR OnQueryDragIcon();
	DECLARE_MESSAGE_MAP()
public:
	afx_msg void OnLbnSelchangeList1();
	afx_msg void OnBnClickedCheck2();
	CString mass_text;
	CString inhomog_text;
	CString inhFWHM_text;
	CString beam_d_text;
	CString crossection_text;
	CString homogRT_text;
	CString homogfinalT_text;
	CString RT_text;
	CString iniT_text;
	CString finalT_text;
	CString cool_t_text;
	CString burn_freq_text;
	CString burnP_text;
	CString burn_t_text;
	CString scanN_text;
	CString t_perscan_text;
	CString min_scan_text;
	CString max_scan_text;
	CString scan_step_text;
	CString recovery_t_text;
	CString recovery_points_text;
	CString cycling_maxT_text;
	CString Ncomplex_text;
	CString Nwells_text;
	CString curvature_text;
	CString barrier_mean_text;
	CString barrier_STD_text;
	CString S_text;
	CString lifetime_text;
	CString stretch_text;
	CString bottom_mean_text;
	CString bottom_STD_text;
	CString Nhires_text;
	CString filename_text;
	CString dist_text;
	CString T_step_text;
	afx_msg void OnBnClickedOk();
	afx_msg void OnBnClickedCancel();
	CString BT_text;
	CString HGK_step_text;
	afx_msg void OnBnClickedEquilibrium();
	CButton Equilibrium;
	CButton Uniform;
	afx_msg void OnBnClickedUniform();
	CButton generate;
	afx_msg void OnBnClickedgenerate();
	CButton balance;
	afx_msg void OnBnClickedbalance();
	CButton bottom;
	afx_msg void OnBnClickedbottom();
	CButton barrier;
	afx_msg void OnBnClickedbarrier();
	CButton analytical;
	afx_msg void OnBnClickedanalytical();
	CButton scan;
	afx_msg void OnBnClickedscan();
	CButton details;
	afx_msg void OnBnClickeddetails();
	CButton one_step;
	afx_msg void OnBnClickedonestep();
	CButton recovery_burn;
	afx_msg void OnBnClickedrecoveryburn();
	afx_msg void OnBnClickedBurn();
	afx_msg void OnBnClickedSmsscan();
	afx_msg void OnBnClickedrecovery();
	afx_msg void OnBnClickedThermocycling();
	CString well_to_STD_ratio;
	afx_msg void OnEnChangeEdit42();
};
