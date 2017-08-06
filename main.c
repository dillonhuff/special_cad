#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
     
#include <poly/assignment.h>
#include <poly/polynomial.h>
#include <poly/variable_order.h>
#include <poly/variable_db.h>
#include <poly/polynomial_context.h>
#include <poly/upolynomial.h>
#include <poly/poly.h>


#include "cad_tree.h"
#include "polynomial_examples.h"

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

// Returns an array of length ps_len.
pl_list all_discriminants(lp_polynomial_t* const * const ps,
			  const size_t ps_len) {
  pl_list p = poly_ptr_list(ps_len);

  for (size_t i = 0; i < ps_len; i++) {
    pl disc = pl_new(lp_polynomial_get_context(ps[i]));

    discriminant(disc, ps[i]);
    
    p[i] = disc;
  }
  
  return p;
}

pl_list all_pairwise_resultants(size_t* num_resultants,
				lp_polynomial_t * const * const ps,
				const size_t ps_len) {
  *num_resultants = (ps_len*(ps_len - 1)) / 2;

  pl_list resultants = poly_ptr_list(*num_resultants);

  size_t total_res = 0;
  for (size_t i = 0; i < ps_len; i++) {
    pl f = ps[i];

    for (size_t j = i + 1; j < ps_len; j++) {
      pl g = ps[j];

      pl res = pl_new(lp_polynomial_get_context(g));
      lp_polynomial_resultant(res, f, g);
      resultants[total_res] = res;
      total_res++;
    }
  }

  assert(total_res == *num_resultants);

  return resultants;
}

int is_duplicate(pl p, lp_polynomial_t* const * const ps, const size_t ps_len) {
  for (int i = 0; i < ps_len; i++) {
    if (lp_polynomial_eq(p, ps[i])) {
      return 1;
    }
  }
  return 0;
}

pl_list remove_duplicates(size_t* no_duplicates_size_ptr,
			  lp_polynomial_t* const * const ps,
			  const size_t ps_len) {
  pl_list nds = poly_ptr_list(ps_len);

  *no_duplicates_size_ptr = 0;

  for (size_t i = 0; i < ps_len; i++) {

    if (!is_duplicate(ps[i], nds, *no_duplicates_size_ptr)) {

      nds[*no_duplicates_size_ptr] = ps[i];
      *no_duplicates_size_ptr += 1;
      
    }
  }

  nds =
    (pl_list)(realloc(nds, sizeof(pl)*(*no_duplicates_size_ptr)));
  
  return nds;
}

pl_list mccallum_projection(size_t* projection_set_size,
			    lp_polynomial_t * const * const ps,
			    const size_t ps_len) {

  // Add coefficients
  size_t num_coeffs = 0;
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

  free(coeffs);

  // Add discriminants
  pl_list discs = all_discriminants(ps, ps_len);

  non_constant_coeffs =
    (pl_list)(realloc(non_constant_coeffs, sizeof(pl)*(*projection_set_size + ps_len)));

  for (size_t i = 0; i < ps_len; i++) {
    if (!lp_polynomial_is_constant(discs[i])) {
      non_constant_coeffs[*projection_set_size] = discs[i];
      *projection_set_size += 1;
    }
  }

  non_constant_coeffs =
    (pl_list)(realloc(non_constant_coeffs, sizeof(pl)*(*projection_set_size)));
  
  free(discs);

  // Add resultants

  size_t num_resultants = 0;
  pl_list resultants = all_pairwise_resultants(&num_resultants, ps, ps_len);

  printf("num_resultants = %zu\n", num_resultants);
  
  non_constant_coeffs =
    (pl_list)(realloc(non_constant_coeffs, sizeof(pl)*(*projection_set_size + num_resultants)));

  for (size_t i = 0; i < num_resultants; i++) {

    printf("resultant %zu = \n", i);
    print_poly(resultants[i]);
    printf("\n");
    
    if (!lp_polynomial_is_constant(resultants[i])) {
      non_constant_coeffs[*projection_set_size] = resultants[i];
      *projection_set_size += 1;
    }
  }

  non_constant_coeffs =
    (pl_list)(realloc(non_constant_coeffs, sizeof(pl)*(*projection_set_size)));

  free(resultants);

  size_t final_proj_size = 0;
  pl_list no_duplicates =
    remove_duplicates(&final_proj_size, non_constant_coeffs, *projection_set_size);

  *projection_set_size = final_proj_size;

  free(non_constant_coeffs);

  return no_duplicates;
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

  assert(proj1_size == 2);

  printf("McCallum projections\n");
  print_poly_list(mc_proj1, proj1_size);

  size_t proj2_size = 0;
  pl_list mc_proj2 = mccallum_projection(&proj2_size, mc_proj1, proj1_size);

  printf("McCallum projection 2\n");
  print_poly_list(mc_proj2, proj2_size);

  assert(proj2_size == 3);

  // Testing lifting

  // Create assignment
  lp_assignment_t* assignment = lp_assignment_new(var_db);


  // Not needed for first lifting
  /* lp_value_t* v_value = lp_value_new(LP_VALUE_INTEGER, &one); */
  /* lp_assignment_set_value(assignment, v, v_value); */

  size_t num_roots = 0;
  lp_value_t* all_roots =
    all_sorted_roots(&num_roots, assignment, mc_proj2, proj2_size);

  printf("# of roots = %zu\n", num_roots);

  for (size_t i = 0; i < num_roots; i++) {

    assert(all_roots[i].type == LP_VALUE_ALGEBRAIC);
    
    lp_algebraic_number_print(&(all_roots[i].value.a), stdout);
    printf("\n");
  }

  size_t num_test_points = 0;
  lp_value_t* all_test_points =
    test_points(&num_test_points, all_roots, num_roots);

  printf("# of test points = %zu\n", num_test_points);

  for (size_t i = 0; i < num_test_points; i++) {

    //assert(all_test_points[i].type == LP_VALUE_ALGEBRAIC);
    
    //lp_algebraic_number_print(&(all_test_points[i].value.a), stdout);
    lp_value_print(&(all_test_points[i]), stdout);

    printf("\n");
  }

  free(all_roots);

  // Projection polynomials to be lifted

  size_t num_projection_sets = 3;

  projection_set* projection_sets =
    (projection_set*)(malloc(sizeof(projection_set)*3));
  projection_sets[0] = make_projection_set(mc_proj2, proj2_size);
  projection_sets[1] = make_projection_set(mc_proj1, proj1_size);
  projection_sets[2] = make_projection_set(&p, 1);


  // Initial empty assignment
  lp_assignment_t* asg = lp_assignment_new(var_db);

  // Create the root of the CAD tree
  cad_cell root = make_cad_cell(NULL, 0, NULL);

  // Actually call CAD lifting
  lift_polynomials(&root, projection_sets, num_projection_sets, asg);

  printf("Final CAD tree\n");
  print_cad_tree(&root);


  size_t total_num_cells = count_cells(&root);
  printf("total # of cells in tree = %zu\n", total_num_cells);

  assert(total_num_cells == 44);
  free(all_test_points);
  lp_assignment_destruct(asg);
}

void test_mccallum_projection_only_resultants() {
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

  assert(proj1_size == 2);

  printf("McCallum projections\n");
  print_poly_list(mc_proj1, proj1_size);

  size_t resultants_size = 0;
  pl_list resultants =
    all_pairwise_resultants(&resultants_size, mc_proj1, proj1_size);

  printf("McCallum projection 2\n");
  printf("# of resultants = %zu\n", resultants_size);
  print_poly_list(resultants, resultants_size);

  lp_integer_destruct(&one);
  lp_integer_destruct(&two);

  for (size_t i = 0; i < proj1_size; i++) {
    lp_polynomial_destruct(mc_proj1[i]);
  }

  free(mc_proj1);

  /* for (size_t i = 0; i < resultants_size; i++) { */
  /*   lp_polynomial_destruct(resultants[i]); */
  /* } */
  
  /* free(resultants); */
}

void test_mccallum_projection() {
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

  assert(proj1_size == 2);

  printf("McCallum projections\n");
  print_poly_list(mc_proj1, proj1_size);

  size_t proj2_size = 0;
  pl_list mc_proj2 = mccallum_projection(&proj2_size, mc_proj1, proj1_size);

  printf("McCallum projection 2\n");
  print_poly_list(mc_proj2, proj2_size);

  assert(proj2_size == 3);

  lp_integer_destruct(&one);
  lp_integer_destruct(&two);

  for (size_t i = 0; i < proj1_size; i++) {
    lp_polynomial_destruct(mc_proj1[i]);
  }

  free(mc_proj1);

  for (size_t i = 0; i < proj2_size; i++) {
    lp_polynomial_destruct(mc_proj2[i]);
  }
  
  free(mc_proj2);
}

void test_conic_sections() {
  lp_variable_db_t* var_db = lp_variable_db_new();

  // Create variable order
  lp_variable_order_t* var_order = lp_variable_order_new();

  lp_polynomial_context_t* ctx =
    lp_polynomial_context_new(lp_Z, var_db, var_order);
  
  lp_polynomial_t** cs = build_2_conic_sections(ctx, var_db, var_order);

  printf("Conic sections\n");
  lp_polynomial_print(cs[0], stdout);
  printf("\n");
  lp_polynomial_print(cs[1], stdout);
  printf("\n");
				       

  clock_t start, end;
  double cpu_time_used;
     
  start = clock();

  size_t projection_set_size = 0;
  pl_list mc_proj1 =
    mccallum_projection(&projection_set_size, cs, 2);

  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;  
  
  printf("Projection set\n");
  print_poly_list(mc_proj1, projection_set_size);

  printf("---------- Time to compute generalized projection = %f\n", cpu_time_used);

  lp_polynomial_context_detach(ctx);
}

void test_constant_conic_sections() {
  lp_variable_db_t* var_db = lp_variable_db_new();

  lp_variable_t x = lp_variable_db_new_variable(var_db, "x");
  lp_variable_t y = lp_variable_db_new_variable(var_db, "y");
  
  // Create variable order
  lp_variable_order_t* var_order = lp_variable_order_new();

  lp_variable_order_push(var_order, x);  
  lp_variable_order_push(var_order, y);
  
  lp_polynomial_context_t* ctx =
    lp_polynomial_context_new(lp_Z, var_db, var_order);

  printf("Specific polynomials\n");

  lp_integer_t* conic_1_coeffs =
    (lp_integer_t*) malloc(sizeof(lp_integer_t)*6);
  conic_1_coeffs[0] = mk_int(1);
  conic_1_coeffs[1] = mk_int(2);
  conic_1_coeffs[2] = mk_int(3);
  conic_1_coeffs[3] = mk_int(4);
  conic_1_coeffs[4] = mk_int(5);
  conic_1_coeffs[5] = mk_int(6);

  pl c1 =
    build_int_coeff_conic_section(ctx, var_db, var_order, conic_1_coeffs, x, y);

  printf("conic section 1 = ");
  print_poly(c1);
  printf("\n");

  lp_integer_t* conic_2_coeffs =
    (lp_integer_t*) malloc(sizeof(lp_integer_t)*6);
  conic_2_coeffs[0] = mk_int(-3);
  conic_2_coeffs[1] = mk_int(5);
  conic_2_coeffs[2] = mk_int(3);
  conic_2_coeffs[3] = mk_int(1);
  conic_2_coeffs[4] = mk_int(2);
  conic_2_coeffs[5] = mk_int(1);

  pl c2 =
    build_int_coeff_conic_section(ctx, var_db, var_order, conic_2_coeffs, x, y);

  printf("conic section 2 = ");
  print_poly(c2);
  printf("\n");

  lp_polynomial_t** cs =
    (lp_polynomial_t**) malloc(sizeof(pl)*2);
  cs[0] = c1;
  cs[1] = c2;

  clock_t start, end;
  double cpu_time_used;

  start = clock();

  size_t projection_set_size = 0;
  pl_list mc_proj1 =
    mccallum_projection(&projection_set_size, cs, 2);

  end = clock();
  cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;  
  
  printf("Projection set\n");
  print_poly_list(mc_proj1, projection_set_size);

  printf("---------- Time to compute generalized projection = %f\n", cpu_time_used);

  // Projection polynomials to be lifted

  size_t num_projection_sets = 2;

  projection_set* projection_sets =
    (projection_set*)(malloc(sizeof(projection_set)*2));
  projection_sets[0] = make_projection_set(mc_proj1, projection_set_size);
  projection_sets[1] = make_projection_set(cs, 2);

  // Initial empty assignment
  lp_assignment_t* asg = lp_assignment_new(var_db);

  // Create the root of the CAD tree
  cad_cell root = make_cad_cell(NULL, 0, NULL);

  // Actually call CAD lifting
  lift_polynomials(&root, projection_sets, num_projection_sets, asg);

  printf("Final CAD tree\n");
  print_cad_tree(&root);

  lp_assignment_destruct(asg);
  
  lp_polynomial_context_detach(ctx);

  free(conic_1_coeffs);
  free(conic_2_coeffs);

}

lp_value_t* dyadic_rational_normalization_error_roots() {
  lp_value_t* roots =
    (lp_value_t*) malloc(sizeof(lp_value_t)*3);

  lp_rational_t five_halves;
  lp_rational_construct_from_int(&five_halves, 5, 2);
  lp_value_construct(&(roots[0]), LP_VALUE_RATIONAL, &five_halves);

  //root[1] = <3*x^2 + (-1*x) + (-1), (-445/1024, -889/2048)>
  //root[2] = <61*x^2 + 8*x + (-8), (-889/2048, -111/256)>

  lp_upolynomial_t* x2 = lp_upolynomial_construct_power(lp_Z, 2, 3);
  lp_upolynomial_t* mx = lp_upolynomial_construct_power(lp_Z, 1, -1);
  lp_upolynomial_t* mone = lp_upolynomial_construct_power(lp_Z, 0, -1);

  lp_upolynomial_t* poly1 = lp_upolynomial_add(x2, mx);
  lp_upolynomial_t* poly = lp_upolynomial_add(poly1, mone);
  size_t deg = lp_upolynomial_degree(poly);

  lp_upolynomial_print(poly, stdout);
  printf("\n");

  lp_dyadic_rational_t num;
  lp_dyadic_rational_construct_from_int(&num, -445, 10);

  lp_dyadic_rational_t denum;
  lp_dyadic_rational_construct_from_int(&denum, -889, 11);

  lp_dyadic_interval_t it;
  lp_dyadic_interval_construct(&it, &num, 1, &denum, 1);

  lp_dyadic_interval_print(&it, stdout);
  printf("\n");

  lp_algebraic_number_t algnum_1;
  lp_algebraic_number_construct(&algnum_1, poly, &it);
  lp_value_construct(&(roots[1]), LP_VALUE_ALGEBRAIC, &algnum_1);

  
  /* lp_algebraic_number_t algnum_2; */
  /* lp_algebraic_number_construct(&algnum_2, poly, &it); */
  /* lp_value_construct(&(roots[2]), LP_VALUE_ALGEBRAIC, &algnum_2); */
  
  return roots;
}

void test_constant_conic_sections_unlifted() {
  lp_variable_db_t* var_db = lp_variable_db_new();

  lp_variable_t x = lp_variable_db_new_variable(var_db, "x");
  lp_variable_t y = lp_variable_db_new_variable(var_db, "y");
  
  // Create variable order
  lp_variable_order_t* var_order = lp_variable_order_new();

  lp_variable_order_push(var_order, x);  
  lp_variable_order_push(var_order, y);
  
  lp_polynomial_context_t* ctx =
    lp_polynomial_context_new(lp_Z, var_db, var_order);

  printf("Specific polynomials\n");

  lp_integer_t* conic_1_coeffs =
    (lp_integer_t*) malloc(sizeof(lp_integer_t)*6);
  conic_1_coeffs[0] = mk_int(1);
  conic_1_coeffs[1] = mk_int(2);
  conic_1_coeffs[2] = mk_int(3);
  conic_1_coeffs[3] = mk_int(4);
  conic_1_coeffs[4] = mk_int(5);
  conic_1_coeffs[5] = mk_int(6);

  pl c1 =
    build_int_coeff_conic_section(ctx, var_db, var_order, conic_1_coeffs, x, y);

  printf("conic section 1 = ");
  print_poly(c1);
  printf("\n");

  lp_integer_t* conic_2_coeffs =
    (lp_integer_t*) malloc(sizeof(lp_integer_t)*6);
  conic_2_coeffs[0] = mk_int(-3);
  conic_2_coeffs[1] = mk_int(5);
  conic_2_coeffs[2] = mk_int(3);
  conic_2_coeffs[3] = mk_int(1);
  conic_2_coeffs[4] = mk_int(2);
  conic_2_coeffs[5] = mk_int(1);

  pl c2 =
    build_int_coeff_conic_section(ctx, var_db, var_order, conic_2_coeffs, x, y);

  printf("conic section 2 = ");
  print_poly(c2);
  printf("\n");

  lp_polynomial_t** cs =
    (lp_polynomial_t**) malloc(sizeof(pl)*2);
  cs[0] = c1;
  cs[1] = c2;

  size_t projection_set_size = 0;
  pl_list mc_proj1 =
    mccallum_projection(&projection_set_size, cs, 2);

  print_poly_list(mc_proj1, projection_set_size);

  // Projection polynomials to be lifted

  printf("about to compute roots\n");
  fflush(stdout);

  lp_assignment_t* asg = lp_assignment_new(var_db);

  /* size_t num_roots = 0; */
  /* lp_value_t* all_roots = //dyadic_rational_normalization_error_roots(); */
  /*   all_sorted_roots(&num_roots, asg, mc_proj1, projection_set_size); */

  size_t num_roots = 3;
  lp_value_t* all_roots = dyadic_rational_normalization_error_roots();

  //all_roots = build_dyadic_normalization_error_roots();

  printf("# of roots = %zu\n", num_roots);

  printf("root[0] with type %u = ", all_roots[0].type);
  lp_value_print(&(all_roots[0]), stdout);
  printf("\n");

  printf("root[1] with type %u = ", all_roots[1].type);
  lp_value_print(&(all_roots[1]), stdout);
  printf("\n");

  printf("root[2] with type %u = ", all_roots[2].type);
  lp_value_print(&(all_roots[2]), stdout);
  printf("\n");
  
  assert(num_roots > 0);

  size_t num_test_points = 0;
  size_t* num_test_points_ptr = &num_test_points;

  *num_test_points_ptr = 2*num_roots + 1;

  // First iteration

  // Construct midpoint
  lp_value_t current = all_roots[0];
  lp_value_t next = all_roots[0 + 1];

  printf("next value = ");
  lp_value_print(&next, stdout);
  printf("\n");

  printf("checking root normalization of\n");
  printf("current = ");
  lp_value_print(&current, stdout);
  printf("\n");
  printf("next = ");
  lp_value_print(&current, stdout);
  printf("\n");

  printf("checking root %d for normalization before between value call\n", 0);
  if (is_algebraic(all_roots[0])) {
    check_normalized(&(all_roots[0].value.a));
  }

  lp_value_t btwn;
  lp_value_construct_none(&btwn);
  lp_value_get_value_between(&current, 1, &next, 1, &btwn);

  printf("checking root %d for normalization after between value\n", 0);
  if (is_algebraic(all_roots[0])) {
    check_normalized(&(all_roots[0].value.a));
  }

  // Second iteration

  // Construct midpoint
  current = all_roots[1];
  next = all_roots[1 + 1];

  printf("next value = ");
  lp_value_print(&next, stdout);
  printf("\n");

  printf("checking root normalization of\n");
  printf("current = ");
  lp_value_print(&current, stdout);
  printf("\n");
  printf("next = ");
  lp_value_print(&current, stdout);
  printf("\n");

  printf("checking root %d for normalization before between value call\n", 1);
  if (is_algebraic(all_roots[1])) {
    check_normalized(&(all_roots[1].value.a));
  }

  lp_value_t btwn1;
  lp_value_construct_none(&btwn1);
  lp_value_get_value_between(&current, 1, &next, 1, &btwn1);

  printf("checking root %d for normalization after between value\n", 1);
  if (is_algebraic(all_roots[1])) {
    check_normalized(&(all_roots[1].value.a));
  }

  // end second iteration

  printf("DONE\n");
  assert(0);
  
  /* lp_value_t current = all_roots[1]; */
  /* lp_value_t next = all_roots[2]; */

  /* printf("checking root normalization of\n"); */
  /* printf("current = "); */
  /* lp_value_print(&current, stdout); */
  /* printf("\n"); */
  /* printf("next = "); */
  /* lp_value_print(&current, stdout); */
  /* printf("\n"); */

  /* if (is_algebraic(current)) { */
  /*   check_normalized(&(current.value.a)); */
  /* } */

  /* lp_value_t btwn; */
  /* lp_value_construct_none(&btwn); */
  /* lp_value_get_value_between(&current, 1, &next, 1, &btwn); */

  /* if (is_algebraic(current)) { */
  /*   check_normalized(&(current.value.a));; */
  /* } */

  // Start of normal code
  /* printf("Checking roots for normalization\n"); */
  /* for (size_t i = 0; i < num_roots; i++) { */
  /*   if (is_algebraic(all_roots[i])) { */
  /*     check_normalized(&(all_roots[i].value.a)); */
  /*   } */
  /* } */

  /* printf("All roots are normalized\n"); */

  /* size_t num_test_points = 0; */
  /* lp_value_t* all_test_points = */
  /*   test_points(&num_test_points, all_roots, num_roots); */

  /* printf("# of test points = %zu\n", num_test_points); */
  
  lp_polynomial_context_detach(ctx);

  free(conic_1_coeffs);
  free(conic_2_coeffs);

}

void test_algebraic_number_copy() {

  lp_dyadic_rational_t d_one;
  lp_dyadic_rational_construct_from_int(&d_one, 2, 0);

  printf("Dyadic rationals\n");

  lp_dyadic_rational_print(&d_one, stdout);
  printf("\n");

  lp_algebraic_number_t a_one;
  lp_algebraic_number_construct_from_dyadic_rational(&a_one, &d_one);

  printf("Algebraic numbers\n");

  lp_algebraic_number_print(&a_one, stdout);
  printf("\n");
  
  lp_algebraic_number_t a_cpy;
  lp_algebraic_number_construct_copy(&a_cpy, &a_one);

  lp_algebraic_number_print(&a_cpy, stdout);
  printf("\n");
  
  lp_algebraic_number_destruct(&a_one);
  lp_algebraic_number_destruct(&a_cpy);

  lp_dyadic_rational_destruct(&d_one);

}

// Actually I think it is this: <3*x^2 + (-1*x) + (-1), (-445/1024, -889/2048)>
// Next value = <61*x^2 + 8*x + (-8), (-889/2048, -111/256)>
// Non-normalized: <3*x^2 + (-1*x) + (-1), (-11529223823181087261/1024, -1729386645466579267/1024)>


void test_algebraic_number_refinement() {
  // Constructing first value
  lp_upolynomial_t* x2 = lp_upolynomial_construct_power(lp_Z, 2, 3);
  lp_upolynomial_t* mx = lp_upolynomial_construct_power(lp_Z, 1, -1);
  lp_upolynomial_t* mone = lp_upolynomial_construct_power(lp_Z, 0, -1);

  lp_upolynomial_t* poly1 = lp_upolynomial_add(x2, mx);
  lp_upolynomial_t* poly = lp_upolynomial_add(poly1, mone);
  size_t deg = lp_upolynomial_degree(poly);

  lp_upolynomial_print(poly, stdout);
  printf("\n");

  lp_dyadic_rational_t num;
  lp_dyadic_rational_construct_from_int(&num, -445, 10);

  lp_dyadic_rational_t denum;
  lp_dyadic_rational_construct_from_int(&denum, -889, 11);

  lp_dyadic_interval_t it;
  lp_dyadic_interval_construct(&it, &num, 1, &denum, 1);

  lp_dyadic_interval_print(&it, stdout);
  printf("\n");

  lp_algebraic_number_t algnum;
  lp_algebraic_number_construct(&algnum, poly, &it);

  lp_algebraic_number_print(&algnum, stdout);
  printf("\n");

  // Constructing next value
  // Next value = <61*x^2 + 8*x + (-8), (-889/2048, -111/256)>
  lp_upolynomial_t* x2_2 = lp_upolynomial_construct_power(lp_Z, 2, 61);
  lp_upolynomial_t* mx_2 = lp_upolynomial_construct_power(lp_Z, 1, 8);
  lp_upolynomial_t* mone_2 = lp_upolynomial_construct_power(lp_Z, 0, -8);

  lp_upolynomial_t* poly1_2 = lp_upolynomial_add(x2_2, mx_2);
  lp_upolynomial_t* poly_2 = lp_upolynomial_add(poly1_2, mone_2);

  lp_upolynomial_print(poly_2, stdout);
  printf("\n");

  lp_dyadic_rational_t num_2;
  lp_dyadic_rational_construct_from_int(&num_2, -889, 11);

  lp_dyadic_rational_t denum_2;
  lp_dyadic_rational_construct_from_int(&denum_2, -111, 8);

  lp_dyadic_interval_t it_2;
  lp_dyadic_interval_construct(&it_2, &num_2, 1, &denum_2, 1);

  lp_dyadic_interval_print(&it, stdout);
  printf("\n");

  lp_algebraic_number_t algnum_2;
  lp_algebraic_number_construct(&algnum_2, poly_2, &it_2);

  lp_algebraic_number_print(&algnum_2, stdout);
  printf("\n");


  lp_value_t val_1;
  lp_value_construct(&val_1, LP_VALUE_ALGEBRAIC, &algnum);

  lp_value_t val_2;
  lp_value_construct(&val_2, LP_VALUE_ALGEBRAIC, &algnum_2);

  printf("Constructing fresh value");

  lp_value_t val;
  lp_value_construct_none(&val);
  lp_value_get_value_between(&val_1, 1, &val_2, 1, &val);

  printf("val_1 after call\n");
  lp_value_print(&val_1, stdout);
  printf("\n");
  
  /* for (size_t i = 0; i < 20; i++) { */
  /*   lp_algebraic_number_refine(&algnum); */

  /*   lp_algebraic_number_print(&algnum, stdout); */
  /*   printf("\n"); */
  /* } */
  
  /* printf("Degree of poly = %zu\n", deg); */

  /* lp_upolynomial_print(poly, stdout); */
  /* printf("\n"); */

  /* lp_algebraic_number_t* roots = */
  /*   (lp_algebraic_number_t*)(malloc(sizeof(lp_algebraic_number_t)*deg)); */

  
}

int main() {
  //test_algebraic_number_refinement();
  /* test_algebraic_number_copy(); */
  /* test_mccallum_projection_only_resultants(); */
  /* isolate_multivariate_roots(); */
  /* test_all_coefficients(); */
  /* test_all_discriminants(); */
  /* test_conic_sections(); */
  //test_constant_conic_sections();
  test_constant_conic_sections_unlifted();

  /* for (int i = 0; i < 100; i++) { */
    
  /*   test_all_discriminants(); */
  /* } */

}
