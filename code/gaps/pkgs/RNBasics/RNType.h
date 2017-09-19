/* Include file GAPS class type identifier */



/* Initialization functions */

int RNInitClassType();
void RNStopClassType();



/* Type definitions */

typedef int RNClassID;



/* Class definition */

class RNClassType {
    public:
        // Constructor functions
        RNClassType(RNClassID id, const char *name);
        RNClassType(const char *name);
	~RNClassType(void);

	// Property functions
	const RNClassID ID(void) const;
	const char *Name(void) const;

    private:
	RNClassID id;
	char *name;
};

	

/* Public constants */

#define RN_UNKNOWN_CLASS_ID (-1)
#define RN_NULL_CLASS_ID (0)



/* Public variables */

extern const RNClassType RNnull_class_type;



/* Useful macros */

#define RN_CLASS_TYPE_DECLARATIONS(__class) \
    virtual const RNClassID ClassID(void) const; \
    virtual const char *ClassName(void) const; \
    virtual const RNClassType& ClassType(void) const; \
    static const RNClassID CLASS_ID(void); \
    static const char * CLASS_NAME(void); \
    static const RNClassType& CLASS_TYPE(void)

#define RN_CLASS_TYPE_DEFINITIONS(__class) \
    const RNClassType __class ## _class_type ( #__class ); \
    const RNClassID __class::ClassID(void) const { return CLASS_ID(); } \
    const char *__class::ClassName(void) const { return CLASS_NAME(); } \
    const RNClassType& __class::ClassType(void) const { return CLASS_TYPE(); } \
    const RNClassID __class::CLASS_ID(void) { return CLASS_TYPE().ID(); } \
    const char * __class::CLASS_NAME(void) { return CLASS_TYPE().Name(); } \
    const RNClassType& __class::CLASS_TYPE(void) { return __class ## _class_type; }




/* Inline functions */

inline const RNClassID RNClassType::
ID(void) const
{
    // Return class id
    return id;
}



inline const char *RNClassType::
Name(void) const
{
    // Return class name
    return name;
}



