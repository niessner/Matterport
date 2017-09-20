#ifndef __RN_MAP__C__
#define __RN_MAP__C__



// Include files

#include "RNBasics.h"



// Member functions 

template <class KeyType, class ValueType>
RNMap<KeyType, ValueType>::
RNMap(void)
{
    // Allocate map
    m = new std::map<KeyType, ValueType, RNMapComparator<KeyType> >();
}



template <class KeyType, class ValueType>
RNMap<KeyType, ValueType>::
RNMap(int (*compare_function)(KeyType, KeyType))
{
    // Allocate map
    RNMapComparator<KeyType> c(compare_function);
    m = new std::map<KeyType, ValueType, RNMapComparator<KeyType> >( c );
}



template <class KeyType, class ValueType>
RNMap<KeyType, ValueType>::
RNMap(const RNMap<KeyType,ValueType>& src)
{
    // Copy map
    m = new std::map<KeyType, ValueType, RNMapComparator<KeyType> >( *(src.m) );
}



template <class KeyType, class ValueType>
RNMap<KeyType, ValueType>::
~RNMap(void)
{
    // Delete map
    delete m;
}




template <class KeyType, class ValueType>
int RNMap<KeyType, ValueType>::
NEntries(void) const
{
    // Return number of entries in map
    return (int) m->size();
}




template <class KeyType, class ValueType>
RNBoolean RNMap<KeyType, ValueType>::
Find(KeyType key, ValueType *value) const
{
    // Return value associated with key, or NULL if not found
    typename std::map<KeyType, ValueType, RNMapComparator<KeyType> >::iterator it = m->find(key);
    if (it == m->end()) return FALSE;
    if (value) *value = it->second;
    return TRUE;
}



template <class KeyType, class ValueType>
void RNMap<KeyType, ValueType>::
Empty(void)
{
    // Remove all entries
    m->clear();
}



template <class KeyType, class ValueType>
void RNMap<KeyType, ValueType>::
Insert(KeyType key, ValueType value)
{
    // Insert entry into map
    (*m)[key] = value;
}



template <class KeyType, class ValueType>
void RNMap<KeyType, ValueType>::
Replace(KeyType key, ValueType value)
{
    // Replace entry in map
    (*m)[key] = value;
}




template <class KeyType, class ValueType>
void RNMap<KeyType, ValueType>::
Remove(KeyType key)
{
    // Remove entry from map
    typename std::map<KeyType, ValueType, RNMapComparator<KeyType> >::iterator it = m->find(key);
    if (it != m->end()) m->erase(it);
}



template <class KeyType, class ValueType>
RNMap<KeyType, ValueType>& RNMap<KeyType, ValueType>::
operator=(const RNMap<KeyType,ValueType>& src)
{
    // Replace map
    delete m;
    m = new std::map<KeyType, ValueType, RNMapComparator<KeyType> >( *(src.m) );
    return *this;
}



#endif
