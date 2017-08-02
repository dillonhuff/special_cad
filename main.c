#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include <poly/assignment.h>
#include <poly/polynomial.h>
#include <poly/variable_order.h>
#include <poly/variable_db.h>
#include <poly/polynomial_context.h>
#include <poly/upolynomial.h>
#include <poly/poly.h>

typedef lp_polynomial_t* pl;
typedef lp_polynomial_t** pl_list;
typedef lp_integer_t lpint;

pl pl_new(const lp_polynomial_context_t* ctx) {
  return lp_polynomial_new(ctx);
}

void print_poly(const lp_polynomial_t* p) {
  lp_polynomial_print(p, stdout);
}

lp_polynomial_t** poly_ptr_list(const size_t len) {
  return (lp_polynomial_t**)(malloc(sizeof(lp_polynomial_t*)*len));
}

void print_poly_list(lp_polynomial_t * const * const ps,
		     const size_t len) {
  for (size_t i = 0; i < len; i++) {
    lp_polynomial_print(ps[i], stdout);
    printf("\n");
  }
}

void coefficients(lp_polynomial_t** coefficients,
		  const lp_polynomial_t* p) {
  for (size_t k = 0; k <= lp_polynomial_degree(p); k++) {
    coefficients[k] = lp_polynomial_new(lp_polynomial_get_context(p));
    lp_polynomial_get_coefficient(coefficients[k], p, k);
  }

}

size_t total_degree(lp_polynomial_t * const * const p,
		    const size_t num_polys) {
  size_t td = 0;

  for (size_t i = 0; i < num_polys; i++) {
    td += lp_polynomial_degree(p[i]) + 1;
  }

  return td;
}

lp_polynomial_t** all_coefficients(size_t* total_num_coefficients,
				   lp_polynomial_t * const * const p,
				   const size_t num_polys) {
  *total_num_coefficients = total_degree(p, num_polys);

  lp_polynomial_t** coeffs = poly_ptr_list(*total_num_coefficients);

  size_t offset = 0;

  for (size_t i = 0; i < num_polys; i++) {
    coefficients(coeffs + offset, p[i]);
    offset += lp_polynomial_degree(p[i]) + 1;
  }

  return coeffs;

}

void print_coefficients(const lp_polynomial_t* p) {
  lp_polynomial_t** coeffs =
    (lp_polynomial_t**)(malloc(sizeof(lp_polynomial_t*)*lp_polynomial_degree(p)));

  coefficients(coeffs, p);
  
  for (size_t k = 0; k <= lp_polynomial_degree(p); k++) {
    lp_polynomial_t* coeff = coeffs[k];
    lp_polynomial_print(coeff, stdout);
    printf("\n");
  }

  free(coeffs);
}

void isolate_multivariate_roots() {

  lp_variable_db_t* var_db = lp_variable_db_new();

  // Create variables
  lp_variable_t u = lp_variable_db_new_variable(var_db, "u");
  lp_variable_t v = lp_variable_db_new_variable(var_db, "v");

  // Create variable order
  lp_variable_order_t* var_order = lp_variable_order_new();
  lp_variable_order_push(var_order, v);  
  lp_variable_order_push(var_order, u);


  lp_polynomial_context_t* ctx =
    lp_polynomial_context_new(lp_Z, var_db, var_order);
  
  lp_integer_t it;
  lp_integer_construct_from_int(lp_Z, &it, 23);

  lp_integer_t one;
  lp_integer_construct_from_int(lp_Z, &one, 1);

  lp_integer_t twelve;
  lp_integer_construct_from_int(lp_Z, &twelve, 12);
  
  
  // Build up polynomial 23u^2 + u from monomials
  lp_polynomial_t* x2 = lp_polynomial_new(ctx);
  lp_polynomial_construct_simple(x2, ctx, &it, u, 2);

  lp_polynomial_print(x2, stdout);
  printf("\n");
  
  lp_polynomial_t* v2 = lp_polynomial_new(ctx);
  lp_polynomial_construct_simple(v2, ctx, &one, v, 2);

  lp_polynomial_print(v2, stdout);
  printf("\n");
  
  lp_polynomial_t* x = lp_polynomial_new(ctx);
  lp_polynomial_construct_simple(x, ctx, &one, u, 1);

  lp_polynomial_print(x, stdout);
  printf("\n");

  lp_polynomial_t* poly1 = lp_polynomial_new(ctx);
  lp_polynomial_add(poly1, x2, x);

  lp_polynomial_t* poly = lp_polynomial_new(ctx);
  lp_polynomial_mul(poly, poly1, v2);

  lp_polynomial_print(poly, stdout);
  printf("\n");

  // Create resultant of 2 polynomials
  lp_polynomial_t* r = lp_polynomial_new(ctx);
  lp_polynomial_construct_simple(r, ctx, &twelve, u, 3);

  lp_polynomial_t* s = lp_polynomial_new(ctx);
  lp_polynomial_construct_simple(s, ctx, &one, u, 0);

  lp_polynomial_t* sr = lp_polynomial_new(ctx);

  lp_polynomial_add(sr, r, s);

  printf("sr = ");
  lp_polynomial_print(sr, stdout);
  printf("\n");

  print_coefficients(sr);
  print_coefficients(poly);
  
  lp_polynomial_t* res = lp_polynomial_new(ctx);
  lp_polynomial_resultant(res, poly, sr);

  printf("Resultant = ");
  lp_polynomial_print(res, stdout);
  printf("\n");

  printf("top variable of sr = %zu\n", lp_polynomial_top_variable(sr));
  printf("top variable of resultant = %zu\n", lp_polynomial_top_variable(res));

  // Create assignment
  lp_assignment_t* assignment = lp_assignment_new(var_db);
  lp_value_t* v_value = lp_value_new(LP_VALUE_INTEGER, &one);
  lp_assignment_set_value(assignment, v, v_value);

  // Isolate roots
  size_t deg = lp_polynomial_degree(poly);
  lp_value_t* roots = malloc(sizeof(lp_value_t)*deg);

  if (!lp_polynomial_is_univariate_m(poly, assignment)) {
    printf("ERROR: ");
    lp_polynomial_print(poly, stdout);
    printf(" is not univariate under the given assignment\n");
    assert(0);
  }

  size_t roots_size = 0;
  lp_polynomial_roots_isolate(poly, assignment, roots, &roots_size);

  printf("# of roots of poly = %zu\n", roots_size);

  for (size_t i = 0; i < roots_size; i++) {
    printf("MULTIVARIATE ROOT = \n");
    printf("Value type = %u\n", roots[i].type);

    assert(roots[i].type == LP_VALUE_ALGEBRAIC);
    
    lp_algebraic_number_print(&(roots[i].value.a), stdout);
    printf("\n");
  }
  
  // Cleanup
  lp_assignment_delete(assignment);
  lp_polynomial_context_detach(ctx);
  lp_variable_db_detach(var_db);
  lp_variable_order_detach(var_order);

}

void isolate_univariate_roots() {

  lp_upolynomial_t* x2 = lp_upolynomial_construct_power(lp_Z, 2, 23);
  lp_upolynomial_t* x = lp_upolynomial_construct_power(lp_Z, 1, 1);

  lp_upolynomial_t* poly = lp_upolynomial_add(x2, x);
  size_t deg = lp_upolynomial_degree(poly);

  printf("Degree of poly = %zu\n", deg);

  lp_upolynomial_print(poly, stdout);
  printf("\n");

  lp_algebraic_number_t* roots =
    (lp_algebraic_number_t*)(malloc(sizeof(lp_algebraic_number_t)*deg));

  size_t roots_size;
  lp_upolynomial_roots_isolate(poly, roots, &roots_size);

  printf("# of roots of poly = %zu\n", roots_size);

  for (size_t i = 0; i < roots_size; i++) {
    printf("ROOT = ");
    lp_algebraic_number_print(&(roots[i]), stdout);
    printf("\n");
  }

}

lp_integer_t mk_int(size_t i) {
  lp_integer_t it;
  lp_integer_construct_from_int(lp_Z, &it, i);

  return it;
}

void test_all_coefficients() {
  lp_variable_db_t* var_db = lp_variable_db_new();

  // Create variables
  lp_variable_t u = lp_variable_db_new_variable(var_db, "u");
  lp_variable_t v = lp_variable_db_new_variable(var_db, "v");

  // Create variable order
  lp_variable_order_t* var_order = lp_variable_order_new();
  lp_variable_order_push(var_order, v);  
  lp_variable_order_push(var_order, u);


  lp_polynomial_context_t* ctx =
    lp_polynomial_context_new(lp_Z, var_db, var_order);
  
  lp_integer_t it;
  lp_integer_construct_from_int(lp_Z, &it, 23);

  lp_integer_t one;
  lp_integer_construct_from_int(lp_Z, &one, 1);

  lp_integer_t twelve;
  lp_integer_construct_from_int(lp_Z, &twelve, 12);
  
  
  // Build up polynomial 23u^2 + u from monomials
  lp_polynomial_t* x2 = lp_polynomial_new(ctx);
  lp_polynomial_construct_simple(x2, ctx, &it, u, 2);

  //lp_polynomial_print(x2, stdout);
  //printf("\n");
  
  lp_polynomial_t* v2 = lp_polynomial_new(ctx);
  lp_polynomial_construct_simple(v2, ctx, &one, v, 2);

  //lp_polynomial_print(v2, stdout);
  //printf("\n");
  
  lp_polynomial_t* x = lp_polynomial_new(ctx);
  lp_polynomial_construct_simple(x, ctx, &one, u, 1);

  //lp_polynomial_print(x, stdout);
  //printf("\n");

  lp_polynomial_t* poly1 = lp_polynomial_new(ctx);
  lp_polynomial_add(poly1, x2, x);

  lp_polynomial_t* poly = lp_polynomial_new(ctx);
  lp_polynomial_mul(poly, poly1, v2);

  //lp_polynomial_print(poly, stdout);
  //printf("\n");

  // Create resultant of 2 polynomials
  lp_polynomial_t* r = lp_polynomial_new(ctx);
  lp_polynomial_construct_simple(r, ctx, &twelve, u, 3);

  lp_polynomial_t* s = lp_polynomial_new(ctx);
  lp_polynomial_construct_simple(s, ctx, &one, u, 0);

  lp_polynomial_t* sr = lp_polynomial_new(ctx);

  lp_polynomial_add(sr, r, s);

  lp_polynomial_t** poly_ptrs =
    poly_ptr_list(2);
  poly_ptrs[0] = sr;
  poly_ptrs[1] = poly;

  lp_polynomial_print(poly_ptrs[0], stdout);
  printf("\n");

  size_t coeffs_len;
  lp_polynomial_t** coeffs =
    all_coefficients(&coeffs_len, poly_ptrs, 2);

  lp_polynomial_print(sr, stdout);
  printf("\n");
  lp_polynomial_print(poly, stdout);
  printf("\n");
  printf("# coeffs = %zu\n", coeffs_len);

  assert(coeffs_len == 7);

  for (size_t i = 0; i < coeffs_len; i++) {
    lp_polynomial_print(coeffs[i], stdout);
    printf("\n");
  }

  size_t coeffs_len2;
  lp_polynomial_t** coeffs2 =
    all_coefficients(&coeffs_len2, coeffs, coeffs_len);

  printf("Coefficients 2\n");
  for (size_t i = 0; i < coeffs_len2; i++) {
    lp_polynomial_print(coeffs2[i], stdout);
    printf("\n");
  }
  
  free(coeffs);
  free(coeffs2);
  
}

void discriminant(lp_polynomial_t* disc, lp_polynomial_t* p) {
  pl deriv = pl_new(lp_polynomial_get_context(p));
  lp_polynomial_derivative(deriv, p);
  printf("deriv = ");
  lp_polynomial_print(deriv, stdout);
  printf("\n");
  lp_polynomial_resultant(disc, p, deriv);
}

pl_list mccallum_projection(size_t* projection_set_size,
			    lp_polynomial_t* const * const ps,
			    const size_t ps_len) {
  size_t num_coeffs;
  pl_list coeffs = all_coefficients(&num_coeffs, ps, ps_len);

  *projection_set_size = 0;
  pl_list non_constant_coeffs = poly_ptr_list(num_coeffs);
  for (size_t i = 0; i < num_coeffs; i++) {

    if (!lp_polynomial_is_constant(coeffs[i])) {

      non_constant_coeffs[*projection_set_size] = coeffs[i];
      *projection_set_size += 1;
      
    }
  }

  non_constant_coeffs =
    (pl_list)(realloc(non_constant_coeffs, sizeof(pl)*(*projection_set_size)));
  return non_constant_coeffs;
}

void test_all_discriminants() {
  lpint two = mk_int(2);
  lpint one = mk_int(1);

  lp_variable_db_t* var_db = lp_variable_db_new();

  // Create variables
  lp_variable_t x1 = lp_variable_db_new_variable(var_db, "x1");
  lp_variable_t x2 = lp_variable_db_new_variable(var_db, "x2");
  lp_variable_t x3 = lp_variable_db_new_variable(var_db, "x3");

  // Create variable order
  lp_variable_order_t* var_order = lp_variable_order_new();
  lp_variable_order_push(var_order, x1);  
  lp_variable_order_push(var_order, x2);
  lp_variable_order_push(var_order, x3);


  lp_polynomial_context_t* ctx =
    lp_polynomial_context_new(lp_Z, var_db, var_order);
  
  pl x1m2sq = lp_polynomial_new(ctx);
  lp_polynomial_construct_simple(x1m2sq, ctx, &one, x1, 1);

  pl twop = lp_polynomial_new(ctx);
  lp_polynomial_construct_simple(twop, ctx, &two, x1, 0);

  lp_polynomial_sub(x1m2sq, x1m2sq, twop);
  lp_polynomial_mul(x1m2sq, x1m2sq, x1m2sq);

  lp_polynomial_print(x1m2sq, stdout);
  printf("\n");

  pl x2m2sq = lp_polynomial_new(ctx);
  lp_polynomial_construct_simple(x2m2sq, ctx, &one, x2, 1);

  lp_polynomial_sub(x2m2sq, x2m2sq, twop);
  lp_polynomial_mul(x2m2sq, x2m2sq, x2m2sq);

  lp_polynomial_print(x2m2sq, stdout);
  printf("\n");

  pl x3m2sq = lp_polynomial_new(ctx);
  lp_polynomial_construct_simple(x3m2sq, ctx, &one, x3, 1);

  lp_polynomial_sub(x3m2sq, x3m2sq, twop);
  lp_polynomial_mul(x3m2sq, x3m2sq, x3m2sq);

  lp_polynomial_print(x3m2sq, stdout);
  printf("\n");

  pl onep = lp_polynomial_new(ctx);
  lp_polynomial_construct_simple(onep, ctx, &one, x1, 0);
  
  pl p = lp_polynomial_new(ctx);
  lp_polynomial_add(p, p, x1m2sq);
  lp_polynomial_add(p, p, x2m2sq);
  lp_polynomial_add(p, p, x3m2sq);
  lp_polynomial_sub(p, p, onep);

  lp_polynomial_print(p, stdout);
  printf("\n");

  size_t proj1_size = 0;
  pl_list mc_proj1 = mccallum_projection(&proj1_size, &p, 1);

  printf("McCallum projections\n");
  print_poly_list(mc_proj1, proj1_size);
  

  size_t coeffs_size;
  pl_list coeffs = all_coefficients(&coeffs_size, &p, 1);

  /* printf("Coefficients\n"); */
  /* print_poly_list(coeffs, coeffs_size); */

  /* pl disc = pl_new(ctx); */

  /* discriminant(disc, p); */

  /* printf("Discriminant = "); */
  /* print_poly(disc); */
  /* printf("\n"); */

}

int main() {
  test_all_coefficients();
  test_all_discriminants();
  //isolate_multivariate_roots();

}
