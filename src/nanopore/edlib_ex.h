#pragma once

#include "edlib_ex_aux.h"

int
Edlib_align(const char* query,
			const int query_size,
			const char* target,
			const int target_size,
			EdlibAlignData* align_data,
			const double error,
			char* query_align,
			char* target_align,
			int* qend,
			int* tend);