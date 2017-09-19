
#pragma once#pragma once

#include "stdafx.h"

#include <vector>
#include <string>
#include <list>



#define X_GLOBAL_APP_STATE_FIELDS \
	X(std::string, s_dataPath) \
	X(std::string, s_sceneFileList) \
	X(unsigned int, s_keyMatchSkip) \
	X(unsigned int, s_keyPadding) \
	X(unsigned int, s_numThresholds) \
	X(unsigned int, s_maxNumImages) \
	X(unsigned int, s_maxNumSensors) \
	X(std::string, s_featureType) \
	X(std::string, s_outputFile) \
	X(unsigned int, s_keyDetectWidth) \
	X(unsigned int, s_keyDetectHeight)
	//X(float, s_siftMatchThresh)


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