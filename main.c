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

lp_value_t* all_sorted_roots(size_t* num_roots_ptr,
			     lp_assignment_t* const assignment,
			     lp_polynomial_t * const * const ps,
			     const size_t ps_len) {
  *num_roots_ptr = 0;
  if (ps_len == 0) {
    return 0;
  }

  lp_value_t* roots =
    (lp_value_t*)(malloc(sizeof(lp_value_t)*lp_polynomial_degree(ps[0])));

  for (size_t pnum = 0; pnum < ps_len; pnum++) {
    // Isolate roots
    size_t deg = lp_polynomial_degree(ps[pnum]);
    lp_value_t* p_roots = malloc(sizeof(lp_value_t)*deg);

    if (!lp_polynomial_is_univariate_m(ps[pnum], assignment)) {
      printf("ERROR: ");
      lp_polynomial_print(ps[pnum], stdout);
      printf(" is not univariate under the given assignment\n");
      assert(0);
    }

    size_t roots_size = 0;
    lp_polynomial_roots_isolate(ps[pnum], assignment, p_roots, &roots_size);
    printf("Initial roots size = %zu\n", roots_size);

    roots = (lp_value_t*)(realloc(roots, sizeof(lp_value_t)*(lp_polynomial_degree(ps[pnum]) + *num_roots_ptr) ));
    for (size_t i = 0; i < roots_size; i++) {
      roots[*num_roots_ptr] = p_roots[i];
      (*num_roots_ptr)++;
    }

    roots = (lp_value_t*)(realloc(roots, sizeof(lp_value_t)*(*num_roots_ptr)));
  }

  qsort(roots, *num_roots_ptr, sizeof(lp_value_t), lp_value_cmp_void);

  return roots;
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
				lp_polynomial_t* const * const ps,
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
			    lp_polynomial_t* const * const ps,
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

    printf("resultant %zu = ", i);
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

size_t is_algebraic(const lp_value_t val) {
  return val.type == LP_VALUE_ALGEBRAIC;
}

void minus_one(lp_algebraic_number_t* res,
	       const lp_algebraic_number_t* const a) {
  assert(a->f == 0);

  lp_dyadic_rational_t tmp;
  lp_dyadic_rational_construct_from_int(&tmp, 1, 0);

  lp_dyadic_rational_t cp;
  lp_dyadic_rational_construct(&cp);
  lp_algebraic_number_get_dyadic_midpoint(a, &cp);

  lp_dyadic_rational_t m_one;
  lp_dyadic_rational_construct(&m_one);
  lp_dyadic_rational_sub(&m_one, &cp, &tmp);

  lp_algebraic_number_construct_from_dyadic_rational(res, &m_one);
  

  /* lp_algebraic_number_t alg_one; */
  /* lp_algebraic_number_construct_from_dyadic_rational(&alg_one, &tmp); */

  /* lp_algebraic_number_t alg_diff; */
  /* lp_algebraic_number_sub(res, a, &alg_one); */

  lp_dyadic_rational_destruct(&m_one);
  lp_dyadic_rational_destruct(&tmp);
  lp_dyadic_rational_destruct(&cp);

  //lp_algebraic_number_destruct(&alg_one);
  
}

lp_value_t* test_points(size_t* num_test_points_ptr,
			lp_value_t const * const all_roots,
			const size_t num_roots) {
  assert(num_roots > 0);

  *num_test_points_ptr = 2*num_roots + 1;

  lp_value_t* test_points =
    (lp_value_t*)(malloc(sizeof(lp_value_t)*(*num_test_points_ptr)));

  lp_algebraic_number_t neg_inf_pt;
  minus_one(&neg_inf_pt, &(all_roots[0].value.a));

  printf("Minus inf point = ");
  lp_algebraic_number_print(&neg_inf_pt, stdout);
  printf("\n");

  // Construct dyadic rational
  /* lp_dyadic_rational_t* fst_point = */
  /*   (lp_dyadic_rational_t*)(malloc(sizeof(lp_dyadic_rational_t))); */
  /* lp_dyadic_rational_from_integer(fst_point, mk_int()); */

  // Construct the rational interval

  //test_points[0] = all_roots[0];//minus_one(all_roots[0].value.a);
  lp_value_construct(&(test_points[0]), LP_VALUE_ALGEBRAIC, &neg_inf_pt);

  size_t index = 1;
  for (size_t i = 0; i < num_roots - 1; i++) {
    test_points[index] = all_roots[i / 2];
    index++;

    // Construct midpoint
    lp_value_t current = all_roots[i];
    lp_value_t next = all_roots[i + 1];

    assert(is_algebraic(current));
    assert(is_algebraic(next));

    /* lp_dyadic_rational_t* dp = */
    /*   (lp_dyadic_rational_t*)(malloc(sizeof(lp_dyadic_rational_t))); */

    lp_dyadic_rational_t cp;
    lp_dyadic_rational_construct(&cp);
    lp_algebraic_number_get_dyadic_midpoint(&(current.value.a), &cp);

    lp_dyadic_rational_t np;
    lp_dyadic_rational_construct(&np);
    lp_algebraic_number_get_dyadic_midpoint(&(next.value.a), &np);

    lp_dyadic_rational_t sum;
    lp_dyadic_rational_construct(&sum);
    lp_dyadic_rational_add(&sum, &cp, &np);

    lp_dyadic_rational_t* mid =
      (lp_dyadic_rational_t*)(malloc(sizeof(lp_dyadic_rational_t)));
    lp_dyadic_rational_construct(mid);
    lp_dyadic_rational_div_2exp(mid, &sum, 1);

    printf("midpoint = ");
    lp_dyadic_rational_print(mid, stdout);
    printf("\n");

    // NOTE: I'm not sure whether I can just use this as a value or whether it
    // needs to be allocated explicitly to survive in the lp_value_t
    lp_algebraic_number_t* mid_a =
      (lp_algebraic_number_t*)(malloc(sizeof(lp_algebraic_number_t)));
    lp_algebraic_number_construct_from_dyadic_rational(mid_a, mid);

    printf("midpoint = ");
    lp_algebraic_number_print(mid_a, stdout);
    printf("\n");
    
    lp_dyadic_rational_destruct(&sum);
    lp_dyadic_rational_destruct(&cp);
    lp_dyadic_rational_destruct(&np);

    /* lp_algebraic_number_t* midpoint = */
    /*   (lp_algebraic_number_t*)(malloc(sizeof(lp_algebraic_number_t))); */
    /* lp_algebraic_number_construct_copy(midpoint, &(current.value.a)); */

    /* lp_algebraic_number_add(midpoint, &(current.value.a), &(next.value.a)); */

    
    lp_value_construct(&(test_points[index]), LP_VALUE_ALGEBRAIC, mid_a);// = all_roots[0];
    index++;
  }

  test_points[*num_test_points_ptr - 2] = all_roots[num_roots - 1];
  test_points[*num_test_points_ptr - 1] = all_roots[0];

  /* // Insert -inf point */
  /* // TODO: Replace this dummy */
  /* test_points[0] = all_roots[0]; */

  /* size_t index = 1; */
  /* // Insert middle points */
  /* for (size_t i = 0; i < num_roots - 1; i++) { */
  /*   test_points[index] = all_roots[0]; */
  /*   index++; */
  /*   test_points[index] = all_roots[0]; */
  /*   index++; */
  /* } */

  /* printf("index = %zu\n", index); */

  /* test_points[index] = all_roots[0]; */

  /* // Insert +inf point */
  /* test_points[*num_test_points_ptr - 1] = all_roots[0]; */

  return test_points;
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

    assert(all_test_points[i].type == LP_VALUE_ALGEBRAIC);
    
    lp_algebraic_number_print(&(all_test_points[i].value.a), stdout);

    printf("\n");
  }
  
  /* free(all_roots); */
  /* free(all_test_points); */

  /* size_t coeffs_size; */
  /* pl_list coeffs = all_coefficients(&coeffs_size, &p, 1); */

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
