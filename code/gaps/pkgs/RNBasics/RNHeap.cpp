// Source file for the heap class 

#ifndef __RN__HEAP__C__
#define __RN__HEAP__C__



// Include files

#include "RNBasics.h"



template <class PtrType>
RNHeap<PtrType>::
RNHeap(RNScalar (*value_callback)(PtrType, void *), 
       PtrType **(*entry_callback)(PtrType, void *), 
       void *callback_data, int least_first)
  : entries(NULL),
    nentries(0),
    nallocated(0),
    value_offset(-1),
    entry_offset(-1),
    value_callback(value_callback),
    entry_callback(entry_callback),
    callback_data(callback_data),
    least_first(least_first)
{
}



template <class PtrType>
RNHeap<PtrType>::
RNHeap(PtrType base, RNScalar *value_ptr, PtrType **entry_ptr, int least_first)
  : entries(NULL),
    nentries(0),
    nallocated(0),
    value_offset(0),
    entry_offset(-1),
    value_callback(NULL),
    entry_callback(NULL),
    callback_data(NULL),
    least_first(least_first)
{
  // Check arguments
  assert(base);
  assert(value_ptr);

  // Compute offsets to data entries in struct referenced by PtrType
  if (value_ptr) value_offset = (unsigned char *) value_ptr - (unsigned char *) base;
  if (entry_ptr) entry_offset = (unsigned char *) entry_ptr - (unsigned char *) base; 

  // Check value offset
  assert(value_offset >= 0);
}



template <class PtrType>
RNHeap<PtrType>::
RNHeap(int value_offset, int entry_offset, int least_first)
  : entries(NULL),
    nentries(0),
    nallocated(0),
    value_offset(value_offset),
    entry_offset(entry_offset),
    value_callback(NULL),
    entry_callback(NULL),
    callback_data(NULL),
    least_first(least_first)
{
  // Check value offset
  assert(value_offset >= 0);
}



template <class PtrType>
RNHeap<PtrType>::
~RNHeap(void)
{
  // Delete entries
  if (entries) delete [] entries;
}



template <class PtrType>
int RNHeap<PtrType>::
IsEmpty(void) const
{
  // Return number of entries
  return (nentries == 0);
}



template <class PtrType>
int RNHeap<PtrType>::
NEntries(void) const
{
  // Return number of entries
  return nentries;
}



template <class PtrType>
PtrType RNHeap<PtrType>::
Kth(int k) const
{
  // Return kth data element
  assert((0 <= k) && (k < nentries));
  return entries[k];
}



template <class PtrType>
PtrType RNHeap<PtrType>::
operator[](int k) const
{
  // Return kth data element
  assert((0 <= k) && (k < nentries));
  return entries[k];
}



template <class PtrType>
PtrType RNHeap<PtrType>::
Peek(void) const
{
  // Check number of entries	
  if (nentries == 0) return NULL;

  // Return head entry
  return entries[0];
}



template <class PtrType>
void RNHeap<PtrType>::
Empty(void)
{
  // Set nentries to zero -- don't deallocate
  Truncate(0);
}



template <class PtrType>
void RNHeap<PtrType>::
Sort(int n)
{
  // Sort the first n entries
  if (n < 0) n = nentries;
  Sort(n, 0, nentries-1);
}



template <class PtrType>
void RNHeap<PtrType>::
Sort(int n, int left, int right)
{
  // Check left and right
  if (left >= right) return;
  if (left >= n) return;

  // Split array based on pivot
  RNScalar pivot_value = Value((left + right) / 2);
  int left_index = left;
  int right_index = right;
  while (left_index < right_index) {
    while ((Value(left_index) <= pivot_value) && (left_index < right_index)) left_index++;
    if (left_index == right_index) break;
    while ((Value(right_index) >= pivot_value) && (left_index < right_index)) right_index--;
    Swap(left_index, right_index);
  }

  // Sort parts recursively
  int pivot_index = left_index;
  Sort(n, left, pivot_index-1);
  Sort(n, pivot_index+1, right);
}



template <class PtrType>
void RNHeap<PtrType>::
Truncate(int n, RNBoolean sort)
{
  // Check nentries
  if (nentries <= n) return;

  // Sort the entries
  if (n > 0) {
    // Sort recursively
    if (sort) Sort(n, 0, nentries - 1);

    // Update entry back pointers 
    if (entry_offset || entry_callback) {
      for (int i = 0; i < n; i++) {
        if (entry_offset >= 0) *((PtrType **) ((unsigned char *) entries[i] + entry_offset)) = &entries[i];
        if (entry_callback) *((PtrType **) (*entry_callback)(entries[i], callback_data)) = &entries[i];
      }
    }
  }

  // Keep only the n best entries
  nentries = n;
}



template <class PtrType>
void RNHeap<PtrType>::
Push(PtrType entry)
{
  // Allocate space for entry
  if (nentries == nallocated) {
    nallocated = (nallocated == 0) ? 1 : 2 * nallocated;
    PtrType *newentries = new PtrType [nallocated];
    for (int i = 0; i < nentries; i++) {
      newentries[i] = entries[i];
      if (entry_offset >= 0) *((PtrType **) ((unsigned char *) newentries[i] + entry_offset)) = &newentries[i];
      if (entry_callback) *((PtrType **) (*entry_callback)(newentries[i], callback_data)) = &newentries[i];
    }
    if (entries) delete [] entries;
    entries = newentries;
  }

  // Put entry into tail
  entries[nentries] = entry;

  // Update entry backpointer
  if (entry_offset >= 0) *((PtrType **) ((unsigned char *) entries[nentries] + entry_offset)) = &entries[nentries];
  if (entry_callback) *((PtrType **) (*entry_callback)(entries[nentries], callback_data)) = &entries[nentries];

  // Increment number of entries
  nentries++;

  // Bubble tail entry up to its rightful spot
  BubbleUp(nentries-1);
}



template <class PtrType>
PtrType RNHeap<PtrType>::
Pop(void)
{
  // Check number of entries	
  if (nentries == 0) return NULL;

  // Get head entry
  PtrType result = entries[0];

  // Update deleted entry backpointer
  if (entry_offset >= 0) *((PtrType **) ((unsigned char *) entries[0] + entry_offset)) = NULL;
  if (entry_callback) *((PtrType **) (*entry_callback)(entries[0], callback_data)) = NULL;

  // Remove head entry, by copying tail over it
  entries[0] = entries[nentries-1];

  // Update new entry[0] backpointer
  if (entry_offset >= 0) *((PtrType **) ((unsigned char *) entries[0] + entry_offset)) = &entries[0];
  if (entry_callback) *((PtrType **) (*entry_callback)(entries[0], callback_data)) = &entries[0];

  // Decrement number of entries
  nentries--;

  // Bubble the head entry down to its rightful spot
  BubbleDown(0);

  // Return original head entry
  return result;
}



template <class PtrType>
void RNHeap<PtrType>::
Update(PtrType entry)
{
  // Search for entry
  PtrType *entryp = NULL;
  if (entry_offset >= 0) entryp = *((PtrType **) ((unsigned char *) entry + entry_offset));
  if (entry_callback) entryp = *((PtrType **) (*entry_callback)(entry, callback_data));
  else {
    // Find entry in heap
    for (int i = 0; i < nentries; i++) {
      if (entries[i] == entry) {
        entryp = &entries[i];
        break;
      }
    }
  }

  // Move data into place
  if (!entryp) return;
  int index = entryp - entries;
  assert((0 <= index) && (index < nentries));
  BubbleUp(BubbleDown(index));
}



template <class PtrType>
void RNHeap<PtrType>::
Remove(PtrType entry)
{
  // Search for entry
  PtrType *entryp = NULL;
  if (entry_offset >= 0) entryp = *((PtrType **) ((unsigned char *) entry + entry_offset));
  if (entry_callback) entryp = *((PtrType **) (*entry_callback)(entry, callback_data));
  else {
    // Find entry in heap
    for (int i = 0; i < nentries; i++) {
      if (entries[i] == entry) {
        entryp = &entries[i];
        break;
      }
    }
  }

  // Move data into place
  if (!entryp) return;
  int index = entryp - entries;
  assert((0 <= index) && (index < nentries));

  // Update deleted entry backpointer
  if (entry_offset >= 0) *((PtrType **) ((unsigned char *) entries[index] + entry_offset)) = NULL;
  if (entry_callback) *((PtrType **) (*entry_callback)(entries[index], callback_data)) = NULL;

  // Remove entry, by copying tail over it
  entries[index] = entries[nentries-1];

  // Update new entry[0] backpointer
  if (entry_offset >= 0) *((PtrType **) ((unsigned char *) entries[index] + entry_offset)) = &entries[index];
  if (entry_callback) *((PtrType **) (*entry_callback)(entries[index], callback_data)) = &entries[index];

  // Decrement number of entries
  nentries--;

  // Bubble the entry down to its rightful spot
  BubbleDown(index);
}



template <class PtrType>
RNScalar RNHeap<PtrType>::
Value(int i) const
{
  // Return value}
  if (value_offset >= 0) {
    return *((RNScalar *) ((unsigned char *) entries[i] + value_offset));
  }
  else if (value_callback) {
    return (*value_callback)(entries[i], callback_data);
  }
  else {
    RNAbort("No heap value callback or offset\n");
    return 0;
  }
}



template <class PtrType>
int RNHeap<PtrType>::
Compare(int i, int j) const
{
  // Get values
  RNScalar value1 = Value(i);
  RNScalar value2 = Value(j);

  // Compare values
  if (least_first) { 
    if (value1 < value2) return -1;
    else if (value1 > value2) return 1;
    else return 0;
  }
  else {
    if (value1 > value2) return -1;
    else if (value1 < value2) return 1;
    else return 0;
  }
}



template <class PtrType>
void RNHeap<PtrType>::
Swap(int i, int j)
{
  // Swap entries
  PtrType swap = entries[i];
  entries[i] = entries[j];
  entries[j] = swap;
  
  // Update entry backpointers
  if (entry_offset >= 0) {
    *((PtrType **) ((unsigned char *) entries[i] + entry_offset)) = &entries[i];
    *((PtrType **) ((unsigned char *) entries[j] + entry_offset)) = &entries[j];
  }
  if (entry_callback) {
    *((PtrType **) (*entry_callback)(entries[i], callback_data)) = &entries[i];
    *((PtrType **) (*entry_callback)(entries[j], callback_data)) = &entries[j];
  }
}



template <class PtrType>
int RNHeap<PtrType>::
BubbleUp(int index)
{
  // Move data up tree to rightful place
  while (index > 0) {
    int parent_index = ((index - 1) >> 1);

    // Compare values
    if (Compare(index, parent_index) >= 0) break;

    // Swap values
    Swap(index, parent_index);

    // Move up tree
    index = parent_index;
  }

  // Return new index
  return index;
}



template <class PtrType>
int RNHeap<PtrType>::
BubbleDown(int index)
{
  // Move data down tree to rightful place
  while (index < nentries) {
    int child_index;
    int child1_index = (index << 1) + 1;
    int child2_index = child1_index + 1;
    RNScalar value = 0;
    RNScalar child_value = 0;
    RNScalar child1_value = 0;
    RNScalar child2_value =0;

    // Check child indices
    if (child1_index >= nentries) {
      break;
    }
    else if (child2_index >= nentries) {
      child_index = child1_index;
    }
    else {
      // Get children values
      if (value_offset >= 0) {
        child1_value = *((RNScalar *) ((unsigned char *) entries[child1_index] + value_offset));
        child2_value = *((RNScalar *) ((unsigned char *) entries[child2_index] + value_offset));
      }
      else if (value_callback) {
        child1_value = (*value_callback)(entries[child1_index], callback_data);
        child2_value = (*value_callback)(entries[child2_index], callback_data);
      }
      else RNAbort("No heap value callback or offset\n");

      // Determine child index
      child_index = child1_index;
      if (least_first) { if (child2_value < child1_value) child_index = child2_index; }
      else { if (child2_value > child1_value) child_index = child2_index; }
    }

    // Get values
    if (value_offset >= 0) {
      value  = *((RNScalar *) ((unsigned char *) entries[index] + value_offset));
      child_value  = *((RNScalar *) ((unsigned char *) entries[child_index] + value_offset));
    }
    else if (value_callback) {
      value = (*value_callback)(entries[index], callback_data);
      child_value = (*value_callback)(entries[child_index], callback_data);
    }
    else RNAbort("No heap value callback or offset\n");

    // Compare values
    if (least_first) { if (value <= child_value) break; }
    else { if (value >= child_value) break; }

    // Swap values
    Swap(index, child_index);

    // Move down tree
    index = child_index;
  }

  // Return new index
  return index;
}



template <class PtrType>
int RNHeap<PtrType>::
IsValid(void)
{
  // Check array
  assert(nallocated >= 0);
  assert(nentries >= 0);
  assert(nentries <= nallocated);
  assert(entries || ((nentries == 0) && (nallocated == 0)));

  // Check backpointers
  if (entry_offset >= 0) {
    for (int i = 0; i < nentries; i++) {
      assert(*((PtrType **) ((unsigned char *) entries[i] + entry_offset)) == &entries[i]);
    }
  }
  else if (entry_callback) {
    for (int i = 0; i < nentries; i++) {
      assert(*((PtrType **) (*entry_callback)(entries[i], callback_data)) == &entries[i]);
    }
  }

  // Check values
  for (int i = 0; i < nentries; i++) {
    // Get  value
    RNScalar value = 0;
    if (value_offset >= 0) value = *((RNScalar *) ((unsigned char *) entries[i] + value_offset));
    else if (value_callback) value = (*value_callback)(entries[i], callback_data);

    // Check child1 value
    int child1_index = (i << 1) + 1;
    if (child1_index < nentries) {
      RNScalar child1_value = 0;
      if (value_offset >= 0) child1_value = *((RNScalar *) ((unsigned char *) entries[child1_index] + value_offset));
      else if (value_callback) child1_value = (*value_callback)(entries[child1_index], callback_data);
      if (least_first) { assert(value <= child1_value); }
      else { assert(value >= child1_value); }
    }

    // Check child2 value
    int child2_index = child1_index + 1;
    if (child2_index < nentries) {
      RNScalar child2_value = 0;
      if (value_offset >= 0) child2_value = *((RNScalar *) ((unsigned char *) entries[child2_index] + value_offset));
      else if (value_callback) child2_value = (*value_callback)(entries[child2_index], callback_data);
      if (least_first) { assert(value <= child2_value); }
      else { assert(value >= child2_value); }
    }
  }

  // Return success
  return 1;
}



#endif









