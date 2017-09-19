/* Include file for the GAPS array class */

#ifndef __RN__ARRAY__H__
#define __RN__ARRAY__H__



/* Library initialization functions */

int RNInitArray();
void RNStopArray();



/* Array entry definition */

typedef void *RNArrayEntry;



/* Array of (void *) class definition */

class RNVArray {
    public:
        // Constructor functions
        RNVArray(void);
        RNVArray(const RNVArray& array);
	~RNVArray(void);

        // Array property functions/operators
	const RNBoolean IsEmpty(void) const;
	const int NAllocated(void) const;
	const int NEntries(void) const;

        // Entry property functions/operators
	const int EntryIndex(const RNArrayEntry *entry) const;
	void *&EntryContents(RNArrayEntry *entry) const;

        // Data access functions/operators
	void *Head(void) const;
	void *Tail(void) const;
	void *Kth(int k) const;
	void *operator[](int k) const;
	void *&operator[](int k);

        // Entry access functions/operators
	RNArrayEntry *HeadEntry(void) const;
	RNArrayEntry *TailEntry(void) const;
	RNArrayEntry *KthEntry(int k) const;
	RNArrayEntry *PrevEntry(const RNArrayEntry *entry) const;
	RNArrayEntry *NextEntry(const RNArrayEntry *entry) const;
	RNArrayEntry *FindEntry(const void *data) const;

        // Insertion functions/operators
	RNArrayEntry *InsertHead(void *data);
	RNArrayEntry *InsertTail(void *data);
	RNArrayEntry *InsertKth(void *data, int k);
	RNArrayEntry *InsertBefore(void *data, RNArrayEntry *entry);
	RNArrayEntry *InsertAfter(void *data, RNArrayEntry *entry);
	RNArrayEntry *Insert(void *data);

        // Removal functions/operators
	void RemoveHead(void);
	void RemoveTail(void);
	void RemoveKth(int k);
	void RemoveEntry(RNArrayEntry *entry);
	void Remove(const void *data);

        // Manipulation functions/operators
	void Empty(RNBoolean deallocate = FALSE);
	void Truncate(int length);
	void Shift(int delta);
	void Shift(int start, int length, int delta);
	void Reverse(void);
	void Reverse(int start, int length);
	void Append(const RNVArray& array);
	void Sort(int (*compare)(const void *data1, const void *data2));
	void BubbleSort(int (*compare)(void *data1, void *data2, void *appl), void *appl);
	void SwapEntries(RNArrayEntry *entry1, RNArrayEntry *entry2);
	void Swap(int i, int j);
	void Resize(int length);
        RNVArray& operator=(const RNVArray& array);

	// Debug function
        RNBoolean IsValid(void) const;

    protected:
	// Internal functions -- do not use these
	RNArrayEntry *InternalInsert(void *data, int k);
	void InternalRemove(int k);

    private:
	RNArrayEntry *entries;
        int nallocated;
	int nentries;
};



/* Inline functions */

#include "RNArray.I"



/* Template class definition */

template <class PtrType>
class RNArray : public RNVArray {
    public:
        // Constructor functions
        RNArray(void) : RNVArray() {};
        RNArray(const RNArray<PtrType>& src) : RNVArray(src) {};

        // Entry property functions/operators
	PtrType& EntryContents(RNArrayEntry *entry) const
            { return (PtrType&) *entry; };

        // Entry access functions/operators
	RNArrayEntry *FindEntry(PtrType data) const
            { return RNVArray::FindEntry((const void *) data); };

        // Data access functions/operators
	PtrType Head(void) const
            { return (PtrType) RNVArray::Head(); };
	PtrType Tail(void) const
            { return (PtrType) RNVArray::Tail(); };
	PtrType Kth(int k) const
            { return (PtrType) RNVArray::Kth(k); };
	PtrType operator[](int k) const
            { return (PtrType) RNVArray::Kth(k); };
	PtrType& operator[](int k)
            { return (PtrType&) RNVArray::operator[](k); };

        // Insertion functions/operators
	RNArrayEntry *InsertHead(PtrType data)
            { return RNVArray::InsertHead((void *) data); };
	RNArrayEntry *InsertTail(PtrType data)
            { return RNVArray::InsertTail((void *) data); };
	RNArrayEntry *InsertKth(PtrType data, int k)
            { return RNVArray::InsertKth((void *) data, k); };
	RNArrayEntry *InsertBefore(PtrType data, RNArrayEntry *entry)
            { return RNVArray::InsertBefore((void *) data, entry); };
	RNArrayEntry *InsertAfter(PtrType data, RNArrayEntry *entry)
            { return RNVArray::InsertAfter((void *) data, entry); };
	RNArrayEntry *Insert(PtrType data)
            { return RNVArray::Insert((void *) data); };

        // Removal functions/operators
	void Remove(const PtrType data)
            { RNVArray::Remove((const void *) data); };
};



#endif






