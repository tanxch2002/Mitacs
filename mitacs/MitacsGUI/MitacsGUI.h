
// MitacsGUI.h: PROJECT_NAME 应用程序的主头文件
//

#pragma once

#ifndef __AFXWIN_H__
	#error "Include'pch.h' before including this file to generate PCH"
#endif

#include "resource.h"		// 主符号


// CMitacsGUIApp:
// For the implementation of this class, see MitacsGUI.cpp
//

class CMitacsGUIApp : public CWinApp
{
public:
	CMitacsGUIApp();

// Rewrite
public:
	virtual BOOL InitInstance();

// accomplish

	DECLARE_MESSAGE_MAP()
};

extern CMitacsGUIApp theApp;
