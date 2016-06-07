#pragma once

#define TEST

#include <cmath>
#include <cstdio>
#include <cstdlib>

#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <map>
#ifdef _OPENMP
#include <omp.h>
#endif

#include <algorithm>
#if defined(__GNUC__) && defined(_OPENMP)
#include <parallel/algorithm>
#endif
#include <stack>

#include <sys/time.h>

#include "avx.h"
#include "namespace.h"
#include "math.h"
#include "vector.h"
#include "matrix.h"

#include "param.h"
#include "EoS.h"
#include "kernel.h"
#include "class.h"
#include "tree.h"

#include "init/RMI.h"

