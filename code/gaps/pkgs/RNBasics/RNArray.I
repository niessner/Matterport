/* Inline functions for GAPS array class */



inline const RNBoolean RNVArray::
IsEmpty(void) const
{
    // Return whether array is empty
    return (nentries == 0);
}



inline const int RNVArray::
NAllocated(void) const
{
    // Return number of entries have allocated memory for
    return nallocated;
}



inline const int RNVArray::
NEntries(void) const
{
    // Return number of entries in array
    return nentries;
}



inline const int RNVArray::
EntryIndex(const RNArrayEntry *entry) const
{
    // Return index of entry in array
    return (entry - entries);
}



inline void *&RNVArray::
EntryContents(RNArrayEntry *entry) const
{
    // Return pointer to data in entry
    return *entry;
}



inline void *RNVArray::
Kth(int k) const
{
    // Return kth data element
    assert(entries);
    assert(nentries >= 0);
    assert((k >= 0) && (k < nentries));
    return entries[k];
}



inline void *RNVArray::
Head(void) const
{
    // Return head data element
    return Kth(0);
}



inline void *RNVArray::
Tail(void) const
{
    // Return tail data element
    return Kth(nentries-1);
}



inline void *RNVArray::
operator[](int k) const
{
    // Return kth data element
    return Kth(k);
}



inline void *&RNVArray::
operator[](int k)
{
    // Return reference to kth data element
    assert(entries);
    assert(nentries >= 0);
    assert((k >= 0) && (k < nentries));
    return entries[k];
}



inline RNArrayEntry *RNVArray::
KthEntry(int k) const
{
    // Return kth entry
    assert(entries);
    assert(nentries >= 0);
    assert((k >= 0) && (k < nentries));
    return &(entries[k]);
}



inline RNArrayEntry *RNVArray::
HeadEntry(void) const
{
    // Return pointer to head entry
    return KthEntry(0);
}



inline RNArrayEntry *RNVArray::
TailEntry(void) const
{
    // Return pointer to tail entry
    return KthEntry(nentries-1);
}



inline RNArrayEntry *RNVArray::
PrevEntry(const RNArrayEntry *entry) const
{
    // Return pointer to previous entry
    int k = EntryIndex(entry);
    assert((k >= 0) && (k < nentries));
    return KthEntry(k-1);
}



inline RNArrayEntry *RNVArray::
NextEntry(const RNArrayEntry *entry) const
{
    // Return pointer to next entry
    int k = EntryIndex(entry);
    assert((k >= 0) && (k < nentries));
    return KthEntry(k+1);
}



inline RNArrayEntry *RNVArray::
InsertHead(void *data)
{
    // Insert data into array at head
    return InternalInsert(data, 0);
}



inline RNArrayEntry *RNVArray::
InsertTail(void *data)
{
    // Insert data into array at tail
    return InternalInsert(data, nentries);
}



inline RNArrayEntry *RNVArray::
InsertKth(void *data, int k)
{
    // Insert data into array in kth position
    return InternalInsert(data, k);
}



inline RNArrayEntry *RNVArray::
InsertBefore(void *data, RNArrayEntry *entry)
{
    // Insert data into array before entry
    return InternalInsert(data, EntryIndex(entry));
}



inline RNArrayEntry *RNVArray::
InsertAfter(void *data, RNArrayEntry *entry)
{
    // Insert data into array after entry
    return InternalInsert(data, EntryIndex(entry)+1);
}



inline RNArrayEntry *RNVArray::
Insert(void *data)
{
    // Insert data into array at default position (tail)
    return InsertTail(data);
}



inline void RNVArray::
RemoveHead(void)
{
    // Remove head entry from array
    InternalRemove(0);
}



inline void RNVArray::
RemoveTail(void)
{
    // Remove tail entry from array
    InternalRemove(nentries-1);
}



inline void RNVArray::
RemoveKth(int k) 
{
    // Remove kth entry from array
    InternalRemove(k);
}



inline void RNVArray::
RemoveEntry(RNArrayEntry *entry)
{
    // Remove data from array
    InternalRemove(EntryIndex(entry));
}



inline void RNVArray::
Remove(const void *data)
{
    // Remove data from array
    RemoveEntry(FindEntry(data));
}



