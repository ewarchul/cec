#ifndef CEC_H
#define CEC_H

#include "interfaces.h"
#include "utils.h"
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

void cec(char *extdatadir, int cecVersion, int problem, double *input, int col,
         int row, double *output, char *suite);

#endif
