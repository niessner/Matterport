/* Include file for GAPS basics module */

#ifndef __RN__BASICS__H__
#define __RN__BASICS__H__



/* Compatability switch include files */

#include "RNBasics/RNCompat.h"



/* External include files */

#include "RNBasics/RNExtern.h"



/* Base class for GAPS modules */

#include "RNBasics/RNBase.h"



/* Error reporting include files */

#include "RNBasics/RNError.h"



/* Memory management include files */

#include "RNBasics/RNMem.h"



/* File management include files */

#include "RNBasics/RNFile.h"



/* Basic bitflags include files */

#include "RNBasics/RNFlags.h"



/* Class type include files */

#include "RNBasics/RNType.h"



/* Math include files */

#include "RNBasics/RNScalar.h"
#include "RNBasics/RNIntval.h"



/* Dynamic array include files */

#include "RNBasics/RNArray.h"
#include "RNBasics/RNQueue.h"
#include "RNBasics/RNHeap.h"
#include "RNBasics/RNMap.h"



/* Graphics utility include files */

#include "RNBasics/RNGrfx.h"
#include "RNBasics/RNRgb.h"



/* OS utility include files */

#include "RNBasics/RNTime.h"



/* SVD stuff */

#include "RNBasics/RNSvd.h"



/* JSON stuff */

#include "RNBasics/json.h"



/* Initialization functions */

int RNInitBasics(void);
void RNStopBasics(void);



#endif







