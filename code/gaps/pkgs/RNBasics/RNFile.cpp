// Source file for file management utilities


// Include files
#include "RNBasics.h"



////////////////////////////////////////////////////////////////////////
// FILE EXISTANCE FUNCTIONS
////////////////////////////////////////////////////////////////////////

RNBoolean
RNFileExists(const char *filename)
{
  // Return whether or not file exists (and is readable)
  FILE *fp = fopen(filename, "rb");
  if (!fp) return FALSE;
  fclose(fp); 
  return TRUE;
}



////////////////////////////////////////////////////////////////////////
// FILE I/O UTILITY FUNCTIONS
////////////////////////////////////////////////////////////////////////

int 
RNFileSeek(FILE *fp, unsigned long long offset, int whence)
{
#if (RN_OS == RN_WINDOWS)
  if (_fseek64(fp, offset, whence) == 0) return 1;
  else return 0;
#else
  // Linux/unix/cygwin etc.
  if (fseeko(fp, offset, whence) == 0) return 1;
  else return 0;
#endif
}



unsigned long long 
RNFileTell(FILE *fp)
{
#if (RN_OS == RN_WINDOWS)
  return _ftell64(fp);
#else
  // Linux/unix/cygwin etc.
  return ftello(fp);
#endif
}



