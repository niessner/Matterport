
#pragma once#pragma once

#include "stdafx.h"

#include <vector>
#include <string>
#include <list>



#define X_GLOBAL_APP_STATE_FIELDS \
	X(std::string, s_srcPath) \
	X(std::string, s_outPath) \
	X(float, s_matchThresh) \
	X(float, s_responseThresh) \
	X(unsigned int, s_maxNumScenes) \
	X(unsigned int, s_maxNumSensFiles) \
	X(unsigned int, s_maxNumImages) \
	X(unsigned int, s_outWidth) \
	X(unsigned int, s_outHeight) \
	X(float, s_depthFilterSigmaD) \
	X(float, s_depthFilterSigmaR) \
	X(float, s_renderDepthMin) \
	X(float, s_renderDepthMax) \
	X(unsigned int, s_maxNumFramesPerScene) \
	X(unsigned int, s_maxNumMatchesPerScene) \
	X(unsigned int, s_maxNumKeysPerFrame)


#ifndef VAR_NAME
#define VAR_NAME(x) #x
#endif

#define checkSizeArray(a, d)( (((sizeof a)/(sizeof a[0])) >= d))

class GlobalAppState
{
public:

#define X(type, name) type name;
	X_GLOBAL_APP_STATE_FIELDS
#undef X

		//! sets the parameter file and reads
	void readMembers(const ParameterFile& parameterFile) {
		m_ParameterFile = parameterFile;
		readMembers();
	}

	//! reads all the members from the given parameter file (could be called for reloading)
	void readMembers() {
#define X(type, name) \
	if (!m_ParameterFile.readParameter(std::string(#name), name)) {MLIB_WARNING(std::string(#name).append(" ").append("uninitialized"));	name = type();}
		X_GLOBAL_APP_STATE_FIELDS
#undef X
 

		m_bIsInitialized = true;
	}

	void print() const {
#define X(type, name) \
	std::cout << #name " = " << name << std::endl;
		X_GLOBAL_APP_STATE_FIELDS
#undef X
	}

	static GlobalAppState& get() {
		static GlobalAppState s;
		return s;
	}


	//! constructor
	GlobalAppState() {
		m_bIsInitialized = false;
	}

	//! destructor
	~GlobalAppState() {
	}

	static void loadGlobalAppState(const std::string& fileNameDescGlobalApp) {
		if (!util::fileExists(fileNameDescGlobalApp)) {
			throw MLIB_EXCEPTION("cannot find parameter filer " + fileNameDescGlobalApp);
		}

		std::cout << VAR_NAME(fileNameDescGlobalApp) << " = " << fileNameDescGlobalApp << std::endl;
		ParameterFile parameterFileGlobalApp(fileNameDescGlobalApp);
		GlobalAppState::get().readMembers(parameterFileGlobalApp);
		GlobalAppState::get().print();
	}

private:
	bool			m_bIsInitialized;
	ParameterFile	m_ParameterFile;
};


typedef GlobalAppState GAS;