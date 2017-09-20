/* Source file for the GAPS array class */



/* Include files */

#include "RNBasics.h"



/* Private constants */

static const int RN_ARRAY_MIN_ALLOCATED = 1;



/* Public functions */

int 
RNInitArray()
{
    /* Return success */
    return TRUE;
}



void 
RNStopArray()
{
}



RNVArray::
RNVArray(void)
    : entries(NULL),
      nallocated(0),
      nentries(0)
{
}



RNVArray::
RNVArray(const RNVArray& array)
    : entries(NULL),
      nallocated(0),
      nentries(0)
{
    // Copy array
    *this = array;
}



 
RNVArray::
~RNVArray(void)
{
    // Free memory for entries
    if (entries) free(entries);
}




RNArrayEntry *RNVArray::
FindEntry(const void *data) const
{
    // Search for entry matching data
    for (int i = 0; i < nentries; i++) 
        if (entries[i] == data) return &entries[i];

    // Entry was not found
    return NULL;
}




RNArrayEntry *RNVArray::
InternalInsert(void *data, int k)
{
    // Make sure there is enough storage for new entry
    Resize(nentries+1);

    // Increment number of entries
    nentries++;

    // Shift entries up one notch 
    if (k < (nentries-1)) Shift(k, 0, 1);

    // Copy data into kth entry
    entries[k] = data;

    // Return entry
    return &entries[k];
}



void RNVArray::
InternalRemove(int k) 
{
    // Shift entries down one notch
    if (k < (nentries-1)) Shift(k+1, 0, -1);

    // Decrement number of entries
    nentries--;
}



void RNVArray::
Empty(RNBoolean deallocate)
{
    // Remove all entries from array
    Truncate(0);

    // Deallocate memory
    if (deallocate) {
        if (entries) delete entries;
        entries = NULL;
        nallocated = 0;
    }
}



void RNVArray::
Truncate(int length)
{
    // Remove tail entries from array
    if (length < nentries) nentries = length;
}



void RNVArray::
Shift(int delta)
{
    // Shift all entries by delta
    Shift(0, 0, delta);
}



void RNVArray::
Shift(int start, int length, int delta)
{
    /* Compute number of entries to shift */
    if ((delta < 0) && (start < -delta)) start = -delta;
    int nshift = nentries - start;
    if (delta > 0) nshift -= delta;
    if (nshift <= 0) return;
    if ((length > 0) && (length < nshift)) nshift = length;

    /* Shift array entries */
    if (delta < 0) {
        for (int i = start; i < (start + nshift); i++) {
	    entries[i+delta]=entries[i];
	}
    }
    else if (delta > 0) {
        for (int i = (start + nshift - 1); i >= start; i--) {
	    entries[i+delta]=entries[i];
	}
    }
}



void RNVArray::
Reverse(void)
{
    // Reverse order of all entries
    Reverse(0, 0);
}



void RNVArray::
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
        Swap(i, j);
    }
}



void RNVArray::
Append(const RNVArray& array)
{
    // Resize first
    Resize(NEntries() + array.NEntries());

    // Insert entries of array
    for (int i = 0; i < array.NEntries(); i++)
	Insert(array.Kth(i));
}



void RNVArray::
Sort(int (*compare)(const void *data1, const void *data2))
{
    // Use qsort
    qsort(entries, nentries, sizeof(void *), compare);
}



void RNVArray::
BubbleSort(int (*compare)(void *data1, void *data2, void *appl), void *appl)
{
    // Sort vector entries
    for (int i = 0; i < NEntries(); i++) {
	for (int j = i+1; j < NEntries(); j++) {
	    if ((*compare)(entries[j], entries[i], appl) < 0) {
	        Swap(i, j);
	    }
	}
    }
}



void RNVArray::
SwapEntries(RNArrayEntry *entry1, RNArrayEntry *entry2)
{
    // Swap entries
    void *tmp = *entry1;
    *entry1 = *entry2;
    *entry2 = tmp;
}



void RNVArray::
Swap(int i, int j)
{
    // Swap entries
    void *tmp = entries[i];
    entries[i] = entries[j];
    entries[j] = tmp;
}



void RNVArray::
Resize(int length)
{
    // Check if length is valid
    assert(length >= nentries);
    assert(nentries <= nallocated);

    // Check if are growing array
    if (length > nallocated) {
	// Adjust length to be next greater power of 2
	int tmplength = RN_ARRAY_MIN_ALLOCATED;
	while (tmplength < length) tmplength *= 2;
	length = tmplength;

	// Allocate new entries
	RNArrayEntry *newentries = NULL;
	if (length > 0) {
	    newentries = (RNArrayEntry *) malloc(length * sizeof(RNArrayEntry));
	    assert(newentries);
	}

	// Copy old entries into new entries
	if (nentries > 0) {
	    assert(entries);
	    assert(newentries);
	    for (int i = 0; i < nentries; i++)
                newentries[i] = entries[i];
	}
 
	// Zero remaining new entries
	if (nentries < length) {
	    assert(newentries);
	    for (int i = nentries; i < length; i++)
                newentries[i] = NULL;
	}

	// Replace entries
	if (entries) free(entries);
	entries = newentries;

	// Update nallocated
	nallocated = length;
    }
}



RNVArray& RNVArray::
operator=(const RNVArray& array)
{
    // Empty array
    Empty();

    // Copy array of entries
    if (array.nentries > 0) {
	// Allocate memory for entries
	Resize(array.nentries);
	
	// Copy entries from array
	for (int i = 0; i < array.nentries; i++) 
            entries[i] = array.entries[i];

	// Update number of entries
	nentries = array.nentries;
    }

    // Return array
    return *this;
}



RNBoolean RNVArray::
IsValid(void) const
{
    // Check invariants
    assert(nentries >= 0);
    assert(nallocated >= 0);
    assert(nentries <= nallocated);
    if (nallocated == 0) { assert(entries == NULL); }
    else { assert(entries != NULL); }

    // Check entries
    for (int i = 0; i < nentries; i++) {
        assert(entries[i]);
    }

    // Return success
    return TRUE;
}






