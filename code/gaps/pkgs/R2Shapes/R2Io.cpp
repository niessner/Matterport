/* Source file for input/output utility */



/* Include files */

#include "R2Shapes/R2Shapes.h"
#include <ctype.h>



RNArray<R2Span *> *
R2ReadSpans(RNArray<R2Span *>& spans, const char *filename)
{
    // Parse input filename extension
    const char *input_extension = strrchr(filename, '.');
    if (!input_extension) {
	printf("Input file has no extension (e.g., .ug).\n");
	return NULL;
    }

    // Read file of appropriate type
    if (!strncmp(input_extension, ".wal", 4)) {
        if (!R2ReadWiseFile(spans, filename)) {
	    RNFail("Unable to open wise file: %s", filename);
	    return NULL;
	}
    }
    else if (!strncmp(input_extension, ".fig", 4)) {
        if (!R2ReadXFigFile(spans, filename)) {
	    RNFail("Unable to open xfig file: %s", filename);
	    return NULL;
	}
    }
    else {
        RNFail("Unrecognized file extension");
	return NULL;
    }

    // Return success
    return &spans;
}



/********************* FORTUNE WISE PARSER *****************************/

RNArray<R2Span *> *
R2ReadWiseFile(RNArray<R2Span *>& spans, const char *filename)
{
    FILE *fp;
    char buffer[1024];
    char *bufferp;
    int line_count = 0;
    char keyword[64], filetype[64];
    int wall_count;
    double x1, y1, z1, x2, y2, z2;
    double dummy;
    int matid;

    // Open file
    fp = fopen(filename, "r");
    if (!fp) {
	RNFail("Unable to open file %s", filename);
	return NULL;
    }

    // Read header
    while (fgets(buffer, 1024, fp)) {
        // Skip comments
        bufferp = buffer;
        line_count++;
	while (isspace(*bufferp)) bufferp++;
	if (*bufferp == '#') continue;
	if (*bufferp == '\0') continue;

	// Parse header: 'filetype wall'
	if (sscanf(bufferp, "%s%s", keyword, filetype) != 2) {
	    RNFail("Syntax error on line %d in file %s", line_count, filename);
	    return NULL;
	}

	// Check keyword
	if (strcmp(keyword, "filetype")) {
	    RNFail("Unrecognized header keyword: %s", keyword);
	    return NULL;
	}

	// Check filetype
	if (strcmp(filetype, "wall")) {
	    RNFail("Unrecognized file type: %s", filetype);
	    return NULL;
	}

	// Done parsing header
	break;
    }

    // Read body
    while (fgets(buffer, 1024, fp)) {
        // Skip comments
        bufferp = buffer;
        line_count++;
	while (isspace(*bufferp)) bufferp++;
	if (*bufferp == '#') continue;
	if (*bufferp == '\0') continue;

	// Parse one span per line: 'wall' i x1 y1 z1 x2 y2 z2 dummy matid
	if (sscanf(bufferp, "%s%d%lf%lf%lf%lf%lf%lf%lf%d",
	    keyword, &wall_count, &x1, &y1, &z1, &x2, &y2, &z2, &dummy, &matid) != 10) {
	    RNFail("Syntax error on line %d in file %s", line_count, filename);
	    return NULL;
	}

	// Check keyword
	if (strcmp(keyword, "wall")) {
	    RNWarning("Unrecognized keyword: %s", keyword);
	    continue;
	}

	// Create span
	if (z1 != z2) {
	    R2Span *span = new R2Span(R2Point(x1, y1), R2Point(x2, y2));
	    spans.Insert(span);
	}
    }
    
    // Close file
    fclose(fp);

    // Return success
    return &spans;
}
    


/********************* XFIG FILE PARSER *****************************/

RNArray<R2Span *> *
R2ReadXFigFile(RNArray<R2Span *>& spans, const char *filename)
{
    FILE *fp;
    char buffer[1024];
    char *bufferp;
    int header_count = 0;
    int type;

    // Open file
    fp = fopen(filename, "r");
    if (!fp) {
	RNFail("Unable to open file %s", filename);
	return NULL;
    }

    // Read header
    while (fgets(buffer, 1024, fp)) {
        // Skip comments
        bufferp = buffer;
	while (isspace(*bufferp)) bufferp++;
	if (*bufferp == '#') continue;
	if (*bufferp == '\0') continue;

	// Increment number of lines parsed in header
        header_count++;

	// Check if done parsing header
	if (header_count == 8) break;
    }

    // Read body
    while (fgets(buffer, 1024, fp)) {
      // Skip comments
      bufferp = buffer;
      while (isspace(*bufferp)) bufferp++;
      if (*bufferp == '#') continue;
      if (*bufferp == '\0') continue;

      // Parse primitive type
      if (sscanf(bufferp, "%d", &type) != 1) {
        RNFail("Syntax error in primitive type in file %s: %s", filename, bufferp);
        return NULL;
      }

      // Check if polyline
      if (type == 2) {
        int npoints;
        int dummy;
        double fdummy;
         
        // Parse polyline fields
        if (sscanf(bufferp, "%d%d%d%d%d%d%d%d%d%lf%d%d%d%d%d%d", 
                   &type, &dummy, &dummy, &dummy, &dummy, &dummy, &dummy, &dummy, 
                   &dummy, &fdummy, &dummy, &dummy, &dummy, &dummy, &dummy, &npoints) != 16) {
          RNFail("Syntax error in primitive fields in file %s: %s", filename, bufferp);
          return NULL;
        }

        // Parse polyline points
        R2Point p1, p2;
        for (int i = 0; i < npoints; i++) {
          // Parse coordinates
          double x, y;
          if (fscanf(fp, "%lf%lf", &x, &y) != 2) {
            RNFail("Syntax error in polyline point in file %s", filename);
            return NULL;
          }

          // Create point
          p2 = R2Point(x, y);

          // Create span
          if (i > 0) {
            R2Span *span = new R2Span(p1, p2);
            spans.Insert(span);
          }

          // Remember point
          p1 = p2;
        }
      }
    }
    
    // Close file
    fclose(fp);

    // Return success
    return &spans;
}
    


