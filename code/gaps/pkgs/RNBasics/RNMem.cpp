/* Source file for GAPS mem utility  */



/* Include files */

#include "RNBasics.h"



int RNInitMem() 
{
    /* Return OK status */
    return TRUE;
}



void RNStopMem()
{
}



#if FALSE

void *operator new(size_t size)
{
    // Allocate size bytes
    void *data = (void *) malloc(size);
    if (!data) RNAbort("Unable to allocate %d bytes", size);
    return data;
}



void *operator new(size_t size, size_t extra)
{
    // Allocate (size+extra) bytes
    void *data = (void *) malloc(size+extra);
    if (!data) RNAbort("Unable to allocate %d bytes", size+extra);
    return (void *) data;
}



void operator delete(void *data)
{
    // Check arguments
    assert(data);

    // Free allocated bytes
    free(data);
}

#endif



void RNSwap(void *node1, void *node2, void *buffer, int size)
{
    // Check arguments
    assert(node1);
    assert(node2);
    assert(size > 0);

    // Find suitable swap buffer 
    void *swap = buffer;
    char swap_buffer[RN_SWAP_BUFFER_SIZE];
    if (buffer == NULL) {
	swap = (void *) swap_buffer;
	if (size > RN_SWAP_BUFFER_SIZE) {
	    swap = (void *) malloc(size);
            if (!swap) {
		RNFail("Unable to allocate swap buffer");
                return;
	    }
	}
    }

    // Swap two buffers 
    RNCopy(node1, swap, size);
    RNCopy(node2, node1, size);
    RNCopy(swap, node2, size);

    // Free swap buffer 
    if (buffer == NULL) {
	if (size > RN_SWAP_BUFFER_SIZE) {
	    free(swap);
	}
    }
}



void RNCopy(const void *src, void *dst, int size)
{
    // Check arguments
    assert(src);
    assert(dst);
    assert(size > 0);

    // Copy buffer
    // bcopy(src, dst, size);
    memcpy(dst, src, size);
}



void RNZero(void *data, int size)
{
    // Check arguments
    assert(data);
    assert(size > 0);

    // Zero buffer
    // bzero(data, size);
    memset(data, 0, size);
}



int RNCompare(const void *src1, const void *src2, int size)
{
    // Check arguments
    assert(src1);
    assert(src2);
    assert(size > 0);

    // Compare buffers
    // return bcmp(src1, src2, size);
    return memcmp(src1, src2, size);
}




long RNMaxMemoryUsage(void)
{
#   if (RN_OS == RN_WINDOWS)
    RNAbort("Not implemented");
    return -1;
#else
    struct rusage usage;
    if (getrusage(RUSAGE_SELF, &usage) == -1) return -1;
    return usage.ru_maxrss;
#   endif
}
