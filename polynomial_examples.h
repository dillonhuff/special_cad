#pragma once

#include <poly/polynomial.h>

#include "utils.h"

lp_polynomial_t** build_2_conic_sections(const lp_polynomial_context_t* ctx,
					 lp_variable_db_t* var_db,
					 lp_variable_order_t* var_order);

lp_polynomial_t*
build_int_coeff_conic_section(const lp_polynomial_context_t* ctx,
			      lp_variable_db_t* var_db,
			      lp_variable_order_t* var_order,
			      lp_integer_t* coefficients,
			      lp_variable_t x,
			      lp_variable_t y);

pl make_plane_polynomial(const lp_polynomial_context_t* ctx,
			 lp_variable_db_t* var_db,
			 lp_variable_order_t* var_order,
			 const lp_variable_t x,
			 const lp_variable_t y,
			 const lp_variable_t z,
			 const lp_variable_t A,
			 const lp_variable_t B,
			 const lp_variable_t C,
			 const lp_variable_t D);

pl make_ellipsoid_polynomial(const lp_polynomial_context_t* ctx,
			     const lp_variable_t x,
			     const lp_variable_t y,
			     const lp_variable_t z,
			     const lp_variable_t E,
			     const lp_variable_t F,
			     const lp_variable_t G,
			     const lp_variable_t H,
			     const lp_variable_t K,
			     const lp_variable_t L);
