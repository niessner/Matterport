/* Inline functions for GAPS queue class */



inline const RNBoolean RNVQueue::
IsEmpty(void) const
{
    // Return whether queue is empty
    return (nentries == 0);
}



inline const int RNVQueue::
NEntries(void) const
{
    // Return number of entries in queue
    return nentries;
}



inline const int RNVQueue::
NAllocated(void) const
{
    // Return number of entries have allocated memory for
    return entries.NEntries();
}



inline const int RNVQueue::
EntryIndex(const RNQueueEntry *entry) const
{
    // Return index of entry in queue
    return KIndex(entries.EntryIndex(entry));
}



inline void *&RNVQueue::
EntryContents(RNQueueEntry *entry) const
{
    // Return pointer to data in entry
    return *entry;
}



inline RNQueueEntry *RNVQueue::
HeadEntry(void) const
{
    // Return pointer to head entry
    return entries.KthEntry(head);
}



inline RNQueueEntry *RNVQueue::
TailEntry(void) const
{
    // Return pointer to tail entry
    return entries.KthEntry(KOffarray(nentries-1));
}



inline RNQueueEntry *RNVQueue::
KthEntry(int k) const
{
    // Return kth entry
    return entries.KthEntry(KOffarray(k));
}



inline RNQueueEntry *RNVQueue::
PrevEntry(const RNQueueEntry *entry) const
{
    // Return pointer to previous entry
    return (entry == entries.HeadEntry()) ? entries.TailEntry() : entries.PrevEntry(entry);
}



inline RNQueueEntry *RNVQueue::
NextEntry(const RNQueueEntry *entry) const
{
    // Return pointer to next entry
    return (entry == entries.TailEntry()) ? entries.HeadEntry() : entries.NextEntry(entry);
}



inline void *RNVQueue::
Kth(int k) const
{
    // Return kth data element
    assert(nentries >= 0);
    assert((k >= 0) && (k < nentries));
    return EntryContents(KthEntry(k));
}



inline void *RNVQueue::
Head(void) const
{
    // Return head data element
    return Kth(0);
}



inline void *RNVQueue::
Tail(void) const
{
    // Return tail data element
    return Kth(nentries-1);
}



inline void *RNVQueue::
operator[](int k) const
{
    // Return kth data element
    return Kth(k);
}



inline RNQueueEntry *RNVQueue::
InsertHead(void *data)
{
    // Insert data into queue at head
    return InternalInsert(data, 0);
}



inline RNQueueEntry *RNVQueue::
InsertTail(void *data)
{
    // Insert data into queue at tail
    return InternalInsert(data, nentries);
}



inline RNQueueEntry *RNVQueue::
InsertKth(void *data, int k)
{
    // Insert data into queue in kth position
    return InternalInsert(data, k);
}



inline RNQueueEntry *RNVQueue::
InsertBefore(void *data, RNQueueEntry *entry)
{
    // Insert data into queue before entry
    return InternalInsert(data, EntryIndex(entry));
}



inline RNQueueEntry *RNVQueue::
InsertAfter(void *data, RNQueueEntry *entry)
{
    // Insert data into queue after entry
    return InternalInsert(data, EntryIndex(entry)+1);
}



inline RNQueueEntry *RNVQueue::
Insert(void *data)
{
    // Insert data into queue at default position (head)
    return InsertTail(data);
}



inline void RNVQueue::
RemoveHead(void)
{
    // Remove head entry from queue
    InternalRemove(0);
}



inline void RNVQueue::
RemoveTail(void)
{
    // Remove tail entry from queue
    InternalRemove(nentries-1);
}



inline void RNVQueue::
RemoveKth(int k) 
{
    // Remove kth entry
    InternalRemove(k);
}



inline void RNVQueue::
RemoveEntry(RNQueueEntry *entry)
{
    // Remove data from queue
    InternalRemove(EntryIndex(entry));
}



inline void RNVQueue::
Remove(const void *data)
{
    // Remove data from queue
    RemoveEntry(FindEntry(data));
}



inline RNQueueEntry *RNVQueue::
Push(void *data)
{
    // Add data to head of queue
    return InsertTail(data);
}



inline void *RNVQueue::
Peek(void)
{
    // Return data from head of queue
    return (IsEmpty()) ? NULL : EntryContents(HeadEntry());
}



