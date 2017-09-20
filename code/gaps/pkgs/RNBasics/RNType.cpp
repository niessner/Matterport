/* Source file for GAPS class types  */



/* Include files */

#include "RNBasics.h"



/* Public variables */

const RNClassType RNnull_class_type(RN_NULL_CLASS_ID, "");



/* Private variables */

static RNClassID RNnclass_ids = 0;



int RNInitClassType() 
{
    /* Return OK status */
    return TRUE;
}



void RNStopClassType()
{
}



RNClassType::
RNClassType(RNClassID id, const char *name)
{
    // Initialize class type
    this->id = id;
    if (!name) this->name = NULL;
    else {
	this->name = new char[strlen(name) + 1];
	assert(this->name);
	strcpy(this->name, name);
    }
}



RNClassType::
RNClassType(const char *name)
{
    // Initialize class type
    id = ++RNnclass_ids;
    if (!name) this->name = NULL;
    else {
	this->name = new char[strlen(name) + 1];
	assert(this->name);
	strcpy(this->name, name);
    }
}



RNClassType::
~RNClassType(void)
{
    // Free name
    if (name) delete [] name;
}

