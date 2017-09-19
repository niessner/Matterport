/* Source file for the GAPS queue class */



/* Include files */

#include "RNBasics.h"



/* Public functions */

int 
RNInitQueue()
{
  /* Return success */
  return TRUE;
}



void 
RNStopQueue()
{
}



RNVQueue::
RNVQueue(void)
  : nentries(0),
    head(0)
{
}



RNVQueue::
RNVQueue(const RNVQueue& queue)
  : entries(queue.entries),
    nentries(queue.nentries),
    head(queue.head)
{
}



RNQueueEntry *RNVQueue::
FindEntry(const void *data) const
{
  // Search for entry matching data
  for (int i = 0; i < nentries; i++) {
    RNQueueEntry *entry = KthEntry(i);
    if (EntryContents(entry) == data) return entry;
  }

  // Entry was not found
  return NULL;
}




RNQueueEntry *RNVQueue::
InternalInsert(void *data, int k)
{
  // Make sure there is enough storage for new entry
  Resize(nentries + 1);

  // Increment number of entries
  nentries++;

  // Check where entry is being inserted
  if (k == 0) {
    // Insert into head by shifting head back one
    head = (head + entries.NEntries() - 1) % entries.NEntries();
  }
  else if (k < (nentries-1)) {
    // Insert into middle by shifting all subsequent entries up one
    Shift(k, 0, 1);
  }

  // Get Kth entry
  RNQueueEntry *entry = KthEntry(k);

  // Insert data in kth entry 
  EntryContents(entry) = data;

  // Return entry
  return entry;
}



void RNVQueue::
InternalRemove(int k) 
{
  // Check where entry is being removed
  if (k == 0) {
    // Remove from head by shifting head forward one
    head = (head + 1) % entries.NEntries();
  }
  else if (k < (nentries-1)) {
    // Remove from middle by shifting all subsequent entries back one
    Shift(k+1, 0, -1);
  }

  // Decrement number of entries
  nentries--;
}



void *RNVQueue::
Pop(void)
{
  // Check if there are any entries 
  if (IsEmpty()) {
    return NULL;
  }
  else {
    // Remove and return data from head of queue
    void *head = Head();
    RemoveHead();
    return head;
  }
}



void RNVQueue::
Empty(void)
{
  // Remove everything from queue
  nentries = 0;
  head = 0;
}



void RNVQueue::
Shift(int delta)
{
  // Shift all entries by delta
  Shift(0, 0, delta);
}



void RNVQueue::
Shift(int start, int length, int delta)
{
  /* Compute number of entries to shift */
  if ((delta < 0) && (start < -delta)) start = -delta;
  int nshift = nentries - start;
  if (delta > 0) nshift -= delta;
  if (nshift <= 0) return;
  if ((length > 0) && (length < nshift)) nshift = length;

  /* Shift queue entries */
  if (delta < 0) {
    for (int i = start; i < (start + nshift); i++) {
      EntryContents(KthEntry(i+delta)) = EntryContents(KthEntry(i));
    }
  }
  else if (delta > 0) {
    for (int i = (start + nshift - 1); i >= start; i--) {
      EntryContents(KthEntry(i+delta)) = EntryContents(KthEntry(i));
    }
  }
}




void RNVQueue::
Reverse(void)
{
  // Reverse order of all entries
  Reverse(0, 0);
}



void RNVQueue::
Reverse(int start, int length)
{
  /* Compute number of entries to reverse */
  int nreverse = nentries - start;
  if (nreverse <= 0) return;
  if ((length > 0) && (length < nreverse)) nreverse = length;
  if (nreverse <= 0) return;

  // Reverse length at start
  int i, j;
  for (i = start, j = start + nreverse - 1; i < j; i++, j--) {
    RNSwap(KthEntry(i), KthEntry(j), NULL, sizeof(void *));
  }
}



void RNVQueue::
Resize(int length)
{
  // Check if are expanding queue
  // Have to rework this if contract queue
  if (length > NAllocated()) {
    // Save number of entries allocated
    int oldnallocated = NAllocated();

    // Resize array of entries - we're reserving all allocated entries
    assert(entries.NAllocated() == entries.NEntries());
    entries.Resize(length);
    int newnallocated = entries.NAllocated();
    for (int i = 0; i < (newnallocated - oldnallocated); i++)
      entries.Insert(NULL);
    assert(entries.NAllocated() == entries.NEntries());
    
    // Check if need to unwrap some entries when array expands
    int oldnwrapped = head + nentries - oldnallocated;
    if (oldnwrapped > 0) {
      int nwrapped = head + nentries - entries.NEntries();
      int nunwrapped = oldnwrapped - nwrapped;
	    
      /* Shift newly unwrapped entries to end */
      if (nunwrapped > 0) entries.Shift(0, nunwrapped, oldnallocated);

      /* Shift still wrapped entries to start */
      if (nwrapped > 0) entries.Shift(nunwrapped, nwrapped, -nunwrapped);
    }
  }
}




RNVQueue& RNVQueue::
operator=(const RNVQueue& queue)
{
  // Assign queue 
  entries = queue.entries;
  nentries = queue.nentries;
  head = queue.head;

  // Return queue
  return *this;
}



RNBoolean RNVQueue::
IsValid(void) const
{
  // Return success
  return TRUE;
}



