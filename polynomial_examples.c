#include "polynomial_examples.h"

#include <stdlib.h>

#include <poly/variable_db.h>

// A_1, B_1, C_1, D_1, E_1, F_1
lp_polynomial_t** build_2_conic_sections(const lp_polynomial_context_t* ctx,
					 lp_variable_db_t* var_db,
					 lp_variable_order_t* var_order) {

  // Create variables
  lp_variable_t x = lp_variable_db_new_variable(var_db, "x");
  lp_variable_t y = lp_variable_db_new_variable(var_db, "y");

  lp_variable_t A_1 = lp_variable_db_new_variable(var_db, "A_1");
  lp_variable_t B_1 = lp_variable_db_new_variable(var_db, "B_1");
  lp_variable_t C_1 = lp_variable_db_new_variable(var_db, "C_1");
  lp_variable_t D_1 = lp_variable_db_new_variable(var_db, "D_1");
  lp_variable_t E_1 = lp_variable_db_new_variable(var_db, "E_1");
  lp_variable_t F_1 = lp_variable_db_new_variable(var_db, "F_1");

  lp_variable_order_push(var_order, A_1);
  lp_variable_order_push(var_order, B_1);
  lp_variable_order_push(var_order, C_1);
  lp_variable_order_push(var_order, D_1);
  lp_variable_order_push(var_order, E_1);
  lp_variable_order_push(var_order, F_1);
  lp_variable_order_push(var_order, x);  
  lp_variable_order_push(var_order, y);

  lpint one = mk_int(1);

  lp_polynomial_t** sections =
    (lp_polynomial_t**) malloc(sizeof(lp_polynomial_t*) * 2);

  sections[0] = pl_new(ctx);
  

  pl xp = pl_new(ctx);
  lp_polynomial_construct_simple(xp, ctx, &one, x, 1);

  pl yp = pl_new(ctx);
  lp_polynomial_construct_simple(yp, ctx, &one, y, 1);
  
  pl a_coeff = pl_new(ctx);
  lp_polynomial_construct_simple(a_coeff, ctx, &one, A_1, 1);

  lp_polynomial_mul(a_coeff, xp, a_coeff);
  lp_polynomial_mul(a_coeff, xp, a_coeff);

  lp_polynomial_add(sections[0], sections[0], a_coeff);
  
  pl b_coeff = pl_new(ctx);
  lp_polynomial_construct_simple(b_coeff, ctx, &one, B_1, 1);

  lp_polynomial_mul(b_coeff, yp, b_coeff);

  lp_polynomial_mul(b_coeff, xp, b_coeff);

  lp_polynomial_add(sections[0], sections[0], b_coeff);

  pl c_coeff = pl_new(ctx);
  lp_polynomial_construct_simple(c_coeff, ctx, &one, C_1, 1);

  lp_polynomial_mul(c_coeff, yp, c_coeff);
  lp_polynomial_mul(c_coeff, yp, c_coeff);

  lp_polynomial_add(sections[0], sections[0], c_coeff);

  pl d_coeff = pl_new(ctx);
  lp_polynomial_construct_simple(d_coeff, ctx, &one, D_1, 1);

  lp_polynomial_mul(d_coeff, xp, d_coeff);

  lp_polynomial_add(sections[0], sections[0], d_coeff);

  pl e_coeff = pl_new(ctx);
  lp_polynomial_construct_simple(e_coeff, ctx, &one, E_1, 1);

  lp_polynomial_mul(e_coeff, yp, e_coeff);
  
  lp_polynomial_add(sections[0], sections[0], e_coeff);

  pl f_coeff = pl_new(ctx);
  lp_polynomial_construct_simple(f_coeff, ctx, &one, F_1, 1);

  lp_polynomial_add(sections[0], sections[0], f_coeff);
  

  sections[1] = pl_new(ctx);

  lp_polynomial_destruct(a_coeff);
  lp_polynomial_destruct(b_coeff);
  lp_polynomial_destruct(c_coeff);
  lp_polynomial_destruct(d_coeff);
  lp_polynomial_destruct(e_coeff);
  lp_polynomial_destruct(xp);
  lp_polynomial_destruct(yp);

  lp_integer_destruct(&one);

  return sections;
  
}
