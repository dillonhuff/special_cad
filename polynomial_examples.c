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

  lp_variable_t A_2 = lp_variable_db_new_variable(var_db, "A_2");
  lp_variable_t B_2 = lp_variable_db_new_variable(var_db, "B_2");
  lp_variable_t C_2 = lp_variable_db_new_variable(var_db, "C_2");
  lp_variable_t D_2 = lp_variable_db_new_variable(var_db, "D_2");
  lp_variable_t E_2 = lp_variable_db_new_variable(var_db, "E_2");
  lp_variable_t F_2 = lp_variable_db_new_variable(var_db, "F_2");
  
  lp_variable_order_push(var_order, A_1);
  lp_variable_order_push(var_order, B_1);
  lp_variable_order_push(var_order, C_1);
  lp_variable_order_push(var_order, D_1);
  lp_variable_order_push(var_order, E_1);
  lp_variable_order_push(var_order, F_1);

  lp_variable_order_push(var_order, A_2);
  lp_variable_order_push(var_order, B_2);
  lp_variable_order_push(var_order, C_2);
  lp_variable_order_push(var_order, D_2);
  lp_variable_order_push(var_order, E_2);
  lp_variable_order_push(var_order, F_2);

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
  
  pl a1_coeff = pl_new(ctx);
  lp_polynomial_construct_simple(a1_coeff, ctx, &one, A_1, 1);

  lp_polynomial_mul(a1_coeff, xp, a1_coeff);
  lp_polynomial_mul(a1_coeff, xp, a1_coeff);

  lp_polynomial_add(sections[0], sections[0], a1_coeff);
  
  pl b1_coeff = pl_new(ctx);
  lp_polynomial_construct_simple(b1_coeff, ctx, &one, B_1, 1);

  lp_polynomial_mul(b1_coeff, yp, b1_coeff);

  lp_polynomial_mul(b1_coeff, xp, b1_coeff);

  lp_polynomial_add(sections[0], sections[0], b1_coeff);

  pl c1_coeff = pl_new(ctx);
  lp_polynomial_construct_simple(c1_coeff, ctx, &one, C_1, 1);

  lp_polynomial_mul(c1_coeff, yp, c1_coeff);
  lp_polynomial_mul(c1_coeff, yp, c1_coeff);

  lp_polynomial_add(sections[0], sections[0], c1_coeff);

  pl d1_coeff = pl_new(ctx);
  lp_polynomial_construct_simple(d1_coeff, ctx, &one, D_1, 1);

  lp_polynomial_mul(d1_coeff, xp, d1_coeff);

  lp_polynomial_add(sections[0], sections[0], d1_coeff);

  pl e1_coeff = pl_new(ctx);
  lp_polynomial_construct_simple(e1_coeff, ctx, &one, E_1, 1);

  lp_polynomial_mul(e1_coeff, yp, e1_coeff);
  
  lp_polynomial_add(sections[0], sections[0], e1_coeff);

  pl f1_coeff = pl_new(ctx);
  lp_polynomial_construct_simple(f1_coeff, ctx, &one, F_1, 1);

  lp_polynomial_add(sections[0], sections[0], f1_coeff);
  

  sections[1] = pl_new(ctx);

  pl a2_coeff = pl_new(ctx);
  lp_polynomial_construct_simple(a2_coeff, ctx, &one, A_2, 1);

  lp_polynomial_mul(a2_coeff, xp, a2_coeff);
  lp_polynomial_mul(a2_coeff, xp, a2_coeff);

  lp_polynomial_add(sections[1], sections[1], a2_coeff);
  
  pl b2_coeff = pl_new(ctx);
  lp_polynomial_construct_simple(b2_coeff, ctx, &one, B_2, 1);

  lp_polynomial_mul(b2_coeff, yp, b2_coeff);

  lp_polynomial_mul(b2_coeff, xp, b2_coeff);

  lp_polynomial_add(sections[1], sections[1], b2_coeff);

  pl c2_coeff = pl_new(ctx);
  lp_polynomial_construct_simple(c2_coeff, ctx, &one, C_2, 1);

  lp_polynomial_mul(c2_coeff, yp, c2_coeff);
  lp_polynomial_mul(c2_coeff, yp, c2_coeff);

  lp_polynomial_add(sections[1], sections[1], c2_coeff);

  pl d2_coeff = pl_new(ctx);
  lp_polynomial_construct_simple(d2_coeff, ctx, &one, D_2, 1);

  lp_polynomial_mul(d2_coeff, xp, d2_coeff);

  lp_polynomial_add(sections[1], sections[1], d2_coeff);

  pl e2_coeff = pl_new(ctx);
  lp_polynomial_construct_simple(e2_coeff, ctx, &one, E_2, 1);

  lp_polynomial_mul(e2_coeff, yp, e2_coeff);
  
  lp_polynomial_add(sections[1], sections[1], e2_coeff);

  pl f2_coeff = pl_new(ctx);
  lp_polynomial_construct_simple(f2_coeff, ctx, &one, F_2, 1);

  lp_polynomial_add(sections[1], sections[1], f2_coeff);

  lp_polynomial_destruct(a1_coeff);
  lp_polynomial_destruct(b1_coeff);
  lp_polynomial_destruct(c1_coeff);
  lp_polynomial_destruct(d1_coeff);
  lp_polynomial_destruct(e1_coeff);
  lp_polynomial_destruct(f1_coeff);
  
  lp_polynomial_destruct(a2_coeff);
  lp_polynomial_destruct(b2_coeff);
  lp_polynomial_destruct(c2_coeff);
  lp_polynomial_destruct(d2_coeff);
  lp_polynomial_destruct(e2_coeff);
  lp_polynomial_destruct(f2_coeff);

  lp_polynomial_destruct(xp);
  lp_polynomial_destruct(yp);

  lp_integer_destruct(&one);

  return sections;
  
}
