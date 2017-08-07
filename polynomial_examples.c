#include "polynomial_examples.h"

#include <stdlib.h>

#include <poly/variable_db.h>

lp_polynomial_t*
build_int_coeff_conic_section(const lp_polynomial_context_t* ctx,
			      lp_variable_db_t* var_db,
			      lp_variable_order_t* var_order,
			      lp_integer_t* coefficients,
			      lp_variable_t x,
			      lp_variable_t y) {
  lp_integer_t one = mk_int(1);
  
  lp_polynomial_t* cs = pl_new(ctx);

  lp_polynomial_t* xs = pl_new(ctx);
  lp_polynomial_construct_simple(xs, ctx, &one, x, 1);

  lp_polynomial_t* ys = pl_new(ctx);
  lp_polynomial_construct_simple(ys, ctx, &one, y, 1);
  
  lp_polynomial_t* xp = pl_new(ctx);
  lp_polynomial_construct_simple(xp, ctx, &(coefficients[3]), x, 1);

  lp_polynomial_add(cs, cs, xp);

  lp_polynomial_t* yp = pl_new(ctx);
  lp_polynomial_construct_simple(yp, ctx, &(coefficients[4]), y, 1);

  lp_polynomial_add(cs, cs, yp);
  
  lp_polynomial_t* xx = pl_new(ctx);
  lp_polynomial_construct_simple(xx, ctx, &(coefficients[0]), x, 2);

  lp_polynomial_add(cs, xx, cs);

  lp_polynomial_t* yy = pl_new(ctx);
  lp_polynomial_construct_simple(yy, ctx, &(coefficients[2]), y, 2);

  lp_polynomial_add(cs, yy, cs);

  lp_polynomial_t* xy = pl_new(ctx);
  lp_polynomial_construct_simple(xy, ctx, &(coefficients[1]), x, 1);
  lp_polynomial_mul(xy, xy, ys);

  lp_polynomial_add(cs, xy, cs);

  lp_polynomial_t* f = pl_new(ctx);
  lp_polynomial_construct_simple(f, ctx, &(coefficients[5]), x, 0);

  lp_polynomial_add(cs, cs, f);

  lp_polynomial_destruct(xs);
  lp_polynomial_destruct(ys);
  lp_polynomial_destruct(xx);
  lp_polynomial_destruct(yy);
  lp_polynomial_destruct(xy);
  lp_polynomial_destruct(xp);
  lp_polynomial_destruct(yp);
  lp_polynomial_destruct(f);


  return cs;
}

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

pl pl_simple_new(const lp_polynomial_context_t* ctx,
		 const size_t a,
		 const lp_variable_t x,
		 const long power) {
  pl p = pl_new(ctx);

  lpint it = mk_int(a);

  lp_polynomial_construct_simple(p, ctx, &it, x, power);

  lp_integer_destruct(&it);
  
  return p;
}

pl make_plane_polynomial(const lp_polynomial_context_t* ctx,
			 lp_variable_db_t* var_db,
			 lp_variable_order_t* var_order,
			 const lp_variable_t x,
			 const lp_variable_t y,
			 const lp_variable_t z,
			 const lp_variable_t A,
			 const lp_variable_t B,
			 const lp_variable_t C,
			 const lp_variable_t D) {
  pl plane = pl_new(ctx);

  
  pl x_term = pl_simple_new(ctx, 1, x, 1);
  pl y_term = pl_simple_new(ctx, 1, y, 1);
  pl z_term = pl_simple_new(ctx, 1, z, 1);
  
  pl a_term = pl_simple_new(ctx, 1, A, 1);

  pl ax = pl_new(ctx);

  lp_polynomial_mul(ax, a_term, x_term);


  lp_polynomial_add(plane, plane, ax);

  pl b_term = pl_simple_new(ctx, 1, B, 1);

  pl bx = pl_new(ctx);

  lp_polynomial_mul(bx, b_term, y_term);

  lp_polynomial_add(plane, plane, bx);

  pl c_term = pl_simple_new(ctx, 1, C, 1);

  pl cx = pl_new(ctx);

  lp_polynomial_mul(cx, c_term, z_term);

  lp_polynomial_add(plane, plane, cx);

  pl d_term = pl_simple_new(ctx, 1, D, 1);

  lp_polynomial_add(plane, plane, d_term);

  pl_delete(ax);
  pl_delete(bx);
  pl_delete(cx);

  pl_delete(a_term);
  pl_delete(b_term);
  pl_delete(c_term);
  pl_delete(d_term);

  pl_delete(x_term);
  pl_delete(y_term);
  pl_delete(z_term);
  
  return plane;
}

pl make_ellipsoid_polynomial(const lp_polynomial_context_t* ctx,
			     const lp_variable_t x,
			     const lp_variable_t y,
			     const lp_variable_t z,
			     const lp_variable_t E,
			     const lp_variable_t F,
			     const lp_variable_t G,
			     const lp_variable_t H,
			     const lp_variable_t K,
			     const lp_variable_t L) {
  pl ellipsoid = pl_new(ctx);

  return ellipsoid;
}
