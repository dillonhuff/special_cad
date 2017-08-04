#pragma once

#include <poly/polynomial.h>

typedef lp_polynomial_t* pl;
typedef lp_polynomial_t** pl_list;
typedef lp_integer_t lpint;

lp_polynomial_t* pl_new(const lp_polynomial_context_t* ctx);
