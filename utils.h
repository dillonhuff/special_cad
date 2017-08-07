#pragma once

#include <poly/polynomial.h>

typedef lp_polynomial_t* pl;
typedef lp_polynomial_t** pl_list;
typedef lp_integer_t lpint;

lp_polynomial_t* pl_new(const lp_polynomial_context_t* ctx);
void pl_delete(pl p);

lp_integer_t mk_int(size_t i);

void print_poly(const lp_polynomial_t* p);

void free_list(void** list, const size_t len);

lp_polynomial_t** poly_ptr_list(const size_t len);

void print_poly_list(lp_polynomial_t * const * const ps,
		     const size_t len);

pl pl_simple_new(const lp_polynomial_context_t* ctx,
		 const size_t a,
		 const lp_variable_t x,
		 const long power);

void set_integer_value(lp_assignment_t* asg, const lp_variable_t x, const int i);
