/* Include file for the GAPS flags class */



/* Initialization functions */

int RNInitFlags();
void RNStopFlags();



/* Class definition */

class RNFlags /* : public RNBase */ {
    public:
        // Constructor functions
        RNFlags(void);
        RNFlags(unsigned long flags);

	// Type conversions
	operator unsigned long(void) const;
	
        // Relationship functions/operators
	int Intersects(const RNFlags flags) const;
	int Contains(const RNFlags flags) const;
  	int operator[](const RNFlags flags) const;

        // Manipulation functions/operators
        void XOR(const RNFlags flags);
        void Add(const RNFlags flags);
        void Remove(const RNFlags flags);
	void Intersect(const RNFlags flags);
        void Reset(unsigned long flags);
  
    private:
	unsigned long flags;
};



/* Flag mask definitions */

#define RN_NO_FLAGS   0x00000000
#define RN_NULL_FLAGS 0x00000000
#define RN_ALL_FLAGS  0xFFFFFFFF



/* Inline functions */

inline RNFlags::
RNFlags(unsigned long flags)
    : flags(flags)
{
}



inline RNFlags::operator 
unsigned long(void) const
{
    // Convert RNFlags to unsigned long
    return flags;
}



inline int RNFlags::
Intersects(const RNFlags flags) const
{
    // Return whether has a property 
    return (this->flags & flags) ? 1 : 0;
}



inline int RNFlags::
Contains(const RNFlags flags) const
{
    // Return whether contains all properties
    return ((this->flags & flags) == flags);
}



inline int RNFlags::
operator[] (const RNFlags flags) const
{
    // Return whether flags intersect
    return this->Intersects(flags);
}



inline void RNFlags::
XOR(const RNFlags flags)
{
    // Union this set of flags with ones passed in
    this->flags ^= flags;
}



inline void RNFlags::
Add(const RNFlags flags)
{
    // Union this set of flags with ones passed in
    this->flags |= flags;
}



inline void RNFlags::
Remove(const RNFlags flags)
{
    // Diff this set of flags with ones passed in
    this->flags &= ~flags;
}



inline void RNFlags::
Intersect(const RNFlags flags)
{
    // Intersect this set of flags with ones passed in
    this->flags &= flags;
}



inline void RNFlags::
Reset(unsigned long flags)
{
  // Set flags
  this->flags = flags;
}




