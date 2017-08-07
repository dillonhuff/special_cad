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
			 lp_variable_order_t* var_order);
