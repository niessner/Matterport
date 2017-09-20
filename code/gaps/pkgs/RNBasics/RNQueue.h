/* Include file for the GAPS queue class */

#ifndef __RN__QUEUE__H__
#define __RN__QUEUE__H__



/* Library initialization functions */

int RNInitQueue();
void RNStopQueue();



/* Entry *definition */

typedef RNArrayEntry RNQueueEntry;



/* Class definition */

class RNVQueue {
    public:
        // Constructor functions
        RNVQueue(void);
        RNVQueue(const RNVQueue& queue);

        // Queue property functions/operators
	const RNBoolean IsEmpty(void) const;
	const int NAllocated(void) const;
	const int NEntries(void) const;

        // Entry *property functions/operators
	const int EntryIndex(const RNQueueEntry *entry) const;
	void *&EntryContents(RNQueueEntry *entry) const;

        // Data access functions/operators
	void *Head(void) const;
	void *Tail(void) const;
	void *Kth(int k) const;
	void *operator[](int k) const;

        // Entry *access functions/operators
	RNQueueEntry *HeadEntry(void) const;
	RNQueueEntry *TailEntry(void) const;
	RNQueueEntry *KthEntry(int k) const;
	RNQueueEntry *PrevEntry(const RNQueueEntry *entry) const;
	RNQueueEntry *NextEntry(const RNQueueEntry *entry) const;
	RNQueueEntry *FindEntry(const void *data) const;

        // Insertion/removal convenience functions
	RNQueueEntry *Push(void *data); 
	void *Pop(void);
	void *Peek(void);

        // Insertion functions/operators
	RNQueueEntry *InsertHead(void *data);
	RNQueueEntry *InsertTail(void *data);
	RNQueueEntry *InsertKth(void *data, int k);
	RNQueueEntry *InsertBefore(void *data, RNQueueEntry *entry);
	RNQueueEntry *InsertAfter(void *data, RNQueueEntry *entry);
	RNQueueEntry *Insert(void *data);

        // Removal functions/operators
	void RemoveHead(void);
	void RemoveTail(void);
	void RemoveKth(int k);
	void RemoveEntry(RNQueueEntry *entry);
	void Remove(const void *data);

        // Manipulation functions/operators
	void Empty(void);
	void Shift(int delta);
	void Shift(int start, int length, int delta);
	void Reverse(void);
	void Reverse(int start, int length);
	void Resize(int length);
        RNVQueue& operator=(const RNVQueue& queue);

	// Debug functions/operators
	RNBoolean IsValid(void) const;

    public:
	// Do not use these
	RNQueueEntry *InternalInsert(void *data, int k);
	void InternalRemove(int k);

    private:
	// Private index re-mapping functions
        int KIndex(int offarray) const { return ((offarray + entries.NAllocated() - head) % entries.NAllocated()); };
        int KOffarray(int index) const { return ((index + head) % entries.NAllocated()); };

    private:
	RNVArray entries;
	int nentries;
	int head;
};



/* Inline functions */

#include "RNQueue.I"



/* Template class definition */

template <class PtrType>
class RNQueue : public RNVQueue {
    public:
        // Constructor functions
        RNQueue(void) : RNVQueue() {};
        RNQueue(const RNQueue<PtrType>& src) : RNVQueue(src) {};

        // Entry property functions/operators
	PtrType& EntryContents(RNQueueEntry *entry) const
            { return (PtrType&) *entry; };

        // Entry access functions/operators
	RNQueueEntry *FindEntry(PtrType data) const
            { return RNVQueue::FindEntry((const void *) data); };

        // Data access functions/operators
	PtrType Head(void) const
            { return (PtrType) RNVQueue::Head(); };
	PtrType Tail(void) const
            { return (PtrType) RNVQueue::Tail(); };
	PtrType Kth(int k) const
            { return (PtrType) RNVQueue::Kth(k); };
	PtrType operator[](int k) const
            { return (PtrType) RNVQueue::Kth(k); };

        // Insertion/removal convenience functions
	void Push(PtrType data)
           { InsertTail(data); };
	PtrType Pop(void)
           { PtrType head = Head(); RemoveHead(); return head; };
	PtrType Peek(void)
           { return Head(); };

        // Insertion functions/operators
	RNQueueEntry *InsertHead(PtrType data)
            { return RNVQueue::InsertHead((void *) data); };
	RNQueueEntry *InsertTail(PtrType data)
            { return RNVQueue::InsertTail((void *) data); };
	RNQueueEntry *InsertKth(PtrType data, int k)
            { return RNVQueue::InsertKth((void *) data, k); };
	RNQueueEntry *InsertBefore(PtrType data, RNQueueEntry *entry)
            { return RNVQueue::InsertBefore((void *) data, entry); };
	RNQueueEntry *InsertAfter(PtrType data, RNQueueEntry *entry)
            { return RNVQueue::InsertAfter((void *) data, entry); };
	RNQueueEntry *Insert(PtrType data)
            { return RNVQueue::Insert((void *) data); };

        // Removal functions/operators
	void Remove(const PtrType data)
            { RNVQueue::Remove((const void *) data); };
};



#endif










