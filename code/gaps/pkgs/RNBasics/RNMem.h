/* Include file for GAPS mem utility */



/* Initialization functions */

int RNInitMem();
void RNStopMem();



/* Memory allocation/deallocation functions */

#if FALSE

void *operator new(size_t size);
void *operator new(size_t size, size_t extra);
void operator delete(void *data);

#endif



/* Standard memory manipulation functions */

#define RN_SWAP_BUFFER_SIZE 1024
void RNSwap(void *node1, void *node2, void *buffer, int size);
void RNCopy(const void *src, void *dst, int size);
void RNZero(void *data, int size);
int RNCompare(const void *src1, const void *src2, int size);



/* Memory usage statistics */

long RNMaxMemoryUsage(void);
