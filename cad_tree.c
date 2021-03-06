#include "cad_tree.h"

#include <assert.h>
#include <stdlib.h>

#include <poly/polynomial.h>
#include <poly/upolynomial.h>

// Copied in from libpoly/src/number/dyadic_rational.h
static inline
int dyadic_rational_is_normalized(const lp_dyadic_rational_t* q) {
  return (mpz_sgn(&q->a) == 0 && q->n == 0) || (mpz_scan1(&q->a, 0) == 0 || q->n == 0);
}

// I think this is the origin of the non-normalized: <3*x^2 + (-1*x) + (-1), (3/4, 1)>
// Actually I think it is this: <3*x^2 + (-1*x) + (-1), (-445/1024, -889/2048)>
// Non-normalized: <3*x^2 + (-1*x) + (-1), (-11529223823181087261/1024, -1729386645466579267/1024)>

void print_dyadic_info(lp_dyadic_rational_t const * q) {
  printf("integer = ");
  lp_integer_print(&(q->a), stdout);
  printf("\n");
  printf("n = %lu\n", q->n);
  printf("is normalized = %d\n", dyadic_rational_is_normalized(q));
  
}

// only a is defined if it is a point
void check_normalized(lp_algebraic_number_t const * a) {
  printf("Checking normalization of ");
  lp_algebraic_number_print(a, stdout);
  printf("\n");
  
  lp_dyadic_interval_t it = a->I;
  printf("Interval = ");
  lp_dyadic_interval_print(&it, stdout);
  printf("\n");

  if (a->f == 0) {
    printf("f == 0\n");
  } else {
    printf("f = ");
    lp_upolynomial_print(a->f, stdout);
    printf("\n");
  }

  printf("a info\n");
  print_dyadic_info(&(it.a));

  if (!(it.is_point)) {
    printf("b info\n");
    print_dyadic_info(&(it.b));
  }

  assert(dyadic_rational_is_normalized(&(it.a)));

  if (!(it.is_point)) {
    assert(dyadic_rational_is_normalized(&(it.b)));
  }

  printf("is normalized\n");
}

void algnum_add(lp_algebraic_number_t* res,
		lp_algebraic_number_t const * a,
		lp_algebraic_number_t const * b) {
  if (a->f && b->f) {
    lp_algebraic_number_add(res, a, b);
    return;
  } else if ((a->f == 0) && (b->f == 0)) {

    lp_dyadic_rational_t tmp_a;
    lp_algebraic_number_get_dyadic_midpoint(a, &tmp_a);

    lp_dyadic_rational_destruct(&tmp_a);

    assert(0);
    return;
  } else {

    printf("ERROR: Mixed rational and algebraic numbers\n");

    printf("a = ");
    //lp_algebraic_number_print(a, stdout);
    printf("\n");

    printf("b = ");
    lp_algebraic_number_print(b, stdout);
    printf("\n");

    assert(0);
  }
}
  //lp_algebraic_number_add(res, a, &tmp_one);

cad_cell make_cad_cell(cad_cell* parent,
		       const size_t num_children,
		       lp_value_t* value) {
  cad_cell c;
  c.parent = parent;
  c.children = (cad_cell*)(malloc(sizeof(cad_cell)*num_children));
  c.num_children = num_children;

  if (value != NULL) {
    c.value = (lp_value_t*) malloc(sizeof(lp_value_t));
    lp_value_construct_none(c.value);
    lp_value_construct_copy(c.value, value);
  } else {
    c.value = NULL;
  }

  return c;
}

projection_set make_projection_set(lp_polynomial_t** polys,
				   size_t length) {
  projection_set s;
  s.polynomials = polys;
  s.length = length;

  return s;
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


void minus_one(lp_algebraic_number_t* res,
	       const lp_algebraic_number_t* const a) {
  if (a->f == 0) {

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
    return;
  }


  assert(a->f != 0);

  printf("minus one a input = ");
  lp_upolynomial_print(a->f, stdout);
  printf("\n");
  
  /* printf("minus one a input = "); */
  /* lp_algebraic_number_print(a, stdout); */
  /* printf("\n"); */

  lp_dyadic_rational_t tmp;
  lp_dyadic_rational_construct_from_int(&tmp, 1, 0);

  lp_dyadic_rational_t tmp2;
  lp_dyadic_rational_construct_from_int(&tmp2, 2, 0);

  lp_dyadic_rational_t tmpm1;
  lp_dyadic_rational_construct_from_int(&tmpm1, -1, 0);
  
  lp_dyadic_interval_t d;
  //lp_dyadic_interval_construct_point(&d, &tmp);
  lp_dyadic_interval_construct(&d, &tmpm1, 1, &tmp2, 1);

  int coeffs[2];
  coeffs[0] = -1;
  coeffs[1] = 1;
  lp_upolynomial_t* f =
    lp_upolynomial_construct_from_int(lp_Z, 1, coeffs);

  printf("Univariate poly = ");
  lp_upolynomial_print(f, stdout);
  printf("\n");
  lp_algebraic_number_t tmp_one;
  lp_algebraic_number_construct(&tmp_one, f, &d);

  printf("Algebraic number 1 = ");
  lp_algebraic_number_print(&tmp_one, stdout);
  printf("\n");

  algnum_add(res, a, &tmp_one);
  //lp_algebraic_number_add(res, a, &tmp_one);

  assert(0);

  lp_dyadic_rational_destruct(&tmp);
  lp_algebraic_number_destruct(&tmp_one);
  lp_dyadic_rational_destruct(&tmp);
  lp_dyadic_rational_destruct(&tmpm1);
  lp_dyadic_rational_destruct(&tmp2);
  lp_dyadic_interval_destruct(&d);
  
}

void plus_one(lp_algebraic_number_t* res,
	       const lp_algebraic_number_t* const a) {
  assert(a->f == 0);

  lp_dyadic_rational_t tmp;
  lp_dyadic_rational_construct_from_int(&tmp, 1, 0);

  lp_dyadic_rational_t cp;
  lp_dyadic_rational_construct(&cp);
  lp_algebraic_number_get_dyadic_midpoint(a, &cp);

  lp_dyadic_rational_t m_one;
  lp_dyadic_rational_construct(&m_one);
  lp_dyadic_rational_add(&m_one, &cp, &tmp);

  lp_algebraic_number_construct_from_dyadic_rational(res, &m_one);

  lp_dyadic_rational_destruct(&m_one);
  lp_dyadic_rational_destruct(&tmp);
  lp_dyadic_rational_destruct(&cp);

}

lp_value_t* test_points(size_t* num_test_points_ptr,
			lp_value_t const * const all_roots,
			const size_t num_roots) {
  if (num_roots == 0) {
    *num_test_points_ptr = 1;

    lp_value_t* test_points =
      (lp_value_t*)(malloc(sizeof(lp_value_t)));

    lp_rational_t tmp;
    lp_rational_construct_from_int(&tmp, 1, 1);

    lp_value_construct(&(test_points[0]), LP_VALUE_RATIONAL, &tmp);

    lp_rational_destruct(&tmp);
    
    return test_points;
  }

  assert(num_roots > 0);

  *num_test_points_ptr = 2*num_roots + 1;

  lp_value_t* test_points =
    (lp_value_t*)(malloc(sizeof(lp_value_t)*(*num_test_points_ptr)));

  lp_value_t neg_inf;
  lp_value_construct(&neg_inf, LP_VALUE_MINUS_INFINITY, NULL);

  lp_value_construct_none(&(test_points[0]));
  lp_value_get_value_between(&neg_inf, 1, &(all_roots[0]), 1, &(test_points[0]));

  lp_value_destruct(&neg_inf);

  printf("Checking first point normalization\n");

  if (is_algebraic(test_points[0])) {
    check_normalized(&(test_points[0].value.a));
  }
  
  size_t index = 1;
  for (size_t i = 0; i < num_roots - 1; i++) {
    test_points[index] = all_roots[i / 2];
    index++;

    // Construct midpoint
    lp_value_t current;// = all_roots[i];
    lp_value_construct_copy(&current, &(all_roots[i]));
    lp_value_t next;// = all_roots[i + 1];
    lp_value_construct_copy(&next, &(all_roots[i + 1]));

    printf("next value = ");
    lp_value_print(&next, stdout);
    printf("\n");

    /* assert(is_algebraic(current)); */
    /* assert(is_algebraic(next)); */

    printf("checking root normalization of\n");
    printf("current = ");
    lp_value_print(&current, stdout);
    printf("\n");
    printf("next = ");
    lp_value_print(&current, stdout);
    printf("\n");

    /* printf("checking root %zu for normalization before between value call\n", i); */
    /* if (is_algebraic(all_roots[i])) { */
    /*   check_normalized(&(all_roots[i].value.a)); */
    /* } */
    
    lp_value_construct_none(&(test_points[index]));
    lp_value_get_value_between(&current, 1, &next, 1, &(test_points[index]));

    lp_value_destruct(&current);
    lp_value_destruct(&next);

    /* printf("checking root %zu for normalization after between value\n", i); */
    /* if (is_algebraic(all_roots[i])) { */
    /*   check_normalized(&(all_roots[i].value.a)); */
    /* } */
    
    /* printf("checking midpoint %zu for normalization\n", index); */
    /* if (is_algebraic(test_points[index])) { */
    /*   check_normalized(&(test_points[index].value.a)); */
    /* } */
    
    index++;
  }

  test_points[*num_test_points_ptr - 2] = all_roots[num_roots - 1];

  lp_value_t pos_inf;
  lp_value_construct(&pos_inf, LP_VALUE_PLUS_INFINITY, NULL);

  lp_value_construct_none(&(test_points[*num_test_points_ptr - 1]));
  lp_value_get_value_between(&(all_roots[num_roots - 1]), 1, &pos_inf, 1, &(test_points[*num_test_points_ptr - 1]));

  lp_value_destruct(&pos_inf);

  /* printf("checking positive endpoint %zu for normalization\n", index); */
  /* if (is_algebraic(test_points[*num_test_points_ptr - 1])) { */
  /*   check_normalized(&(test_points[*num_test_points_ptr - 1].value.a)); */
  /* } */
  
  printf("Testing all test points for normalization\n");

  // Check that all points are normalized
  for (size_t i = 0; i < *num_test_points_ptr; i++) {

    if (is_algebraic(test_points[i])) {
      printf("Checking %zu for normalization\n", i);
      check_normalized(&(test_points[i].value.a));
    }
  }

  printf("All test points are normalized\n");

  return test_points;
}

void lift_polynomials(cad_cell* root,
		      projection_set* projection_sets,
		      const size_t num_projection_sets,
		      lp_assignment_t* asg) {
  printf("LIFTING, # of projection sets left = %zu\n", num_projection_sets);
  printf("CURRENT ASSIGNMENT\n");
  lp_assignment_print(asg, stdout);
  printf("\n");

  if (num_projection_sets == 0) {
    printf("DONE\n");
    return;
  }

  projection_set p = projection_sets[0];

  assert(p.length > 0);
  lp_variable_t next_var = lp_polynomial_top_variable(p.polynomials[0]);

  projection_set* sets_left = projection_sets + 1;

  printf("about to compute roots\n");
  fflush(stdout);

  size_t num_roots = 0;
  lp_value_t* all_roots =
    all_sorted_roots(&num_roots, asg, p.polynomials, p.length);

  printf("# of roots = %zu\n", num_roots);

  printf("Checking roots for normalization\n");
  for (size_t i = 0; i < num_roots; i++) {
    if (is_algebraic(all_roots[i])) {
      check_normalized(&(all_roots[i].value.a));
    }
  }

  printf("All roots are normalized\n");

  size_t num_test_points = 0;
  lp_value_t* all_test_points =
    test_points(&num_test_points, all_roots, num_roots);

  printf("# of test points = %zu\n", num_test_points);

  cad_cell* children =
    (cad_cell*)(malloc(sizeof(cad_cell) * num_test_points));

  fflush(stdout);

  for (size_t i = 0; i < num_test_points; i++) {

    //assert(all_test_points[i].type == LP_VALUE_ALGEBRAIC);
    //lp_algebraic_number_print(&(all_test_points[i].value.a), stdout);
    lp_value_print(&(all_test_points[i]), stdout);
    printf("\n");


    // TODO: Actuall construct new cells
    children[i] = make_cad_cell(root, 0, &(all_test_points[i]));

    /* printf("lp value type = %u\n", all_test_points[i].type); */
    /* printf("lp value = "); */
    /* lp_value_print(&(all_test_points[i]), stdout); */
    /* printf("\n"); */

    /* if (is_algebraic(all_test_points[i])) { */
    /*   printf("IS ALGEBRAIC\n"); */
    /*   lp_dyadic_interval_t it = all_test_points[i].value.a.I; */
    /*   printf("Interval = "); */
    /*   lp_dyadic_interval_print(&it, stdout); */
    /*   printf("\n"); */

    /*   if (all_test_points[i].value.a.f == 0) { */
    /* 	printf("f == 0\n"); */
    /*   } else { */
    /* 	printf("f = "); */
    /* 	lp_upolynomial_print(all_test_points[i].value.a.f, stdout); */
    /* 	printf("\n"); */
    /*   } */

    /*   check_normalized(&(all_test_points[i].value.a)); */

    /* } else { */
    /*   printf("NOT ALGEBRAIC\n"); */
    /* } */

    lp_assignment_set_value(asg, next_var, &(all_test_points[i]));

    lift_polynomials(&(children[i]), sets_left, num_projection_sets - 1, asg);

    lp_assignment_set_value(asg, next_var, NULL);

    //lp_value_destruct(&fresh_value);
  }

  root->children = children;
  root->num_children = num_test_points;

  /* for (size_t i = 0; i < num_test_points; i++) { */
  /*   lp_value_destruct(&(all_test_points[i])); */
  /* } */
  free(all_test_points);

  if (all_roots) {
    free(all_roots);
  }
  
}

void tab_out(const size_t level) {
  for (size_t i = 0; i < level; i++) {
    printf("\t");
  }
}

void print_cad_tree_rec(const size_t level, cad_cell* root) {
  if (root->value == NULL) {
    tab_out(level);
    printf("ROOT\n");
  } else {
    //assert(is_algebraic(*(root->value)));

    tab_out(level);
    printf("CELL VALUE = ");
    //lp_algebraic_number_print(&((root->value)->value.a), stdout);
    lp_value_print(root->value, stdout);
    printf("\n");
  }

  for (size_t i = 0; i < root->num_children; i++) {
    
    print_cad_tree_rec(level + 1, &(root->children[i]));
  }
}

void print_cad_tree(cad_cell* root) {
  print_cad_tree_rec(0, root);
}

size_t is_algebraic(const lp_value_t val) {
  return val.type == LP_VALUE_ALGEBRAIC;
}

size_t is_rational(const lp_value_t val) {
  return val.type == LP_VALUE_RATIONAL;
}

size_t count_cells(cad_cell const * root) {
  size_t ncells = 1;
  for (size_t i = 0; i < root->num_children; i++) {
    ncells += count_cells(&(root->children[i]));
  }
  return ncells;
}

void cad_tree_destruct(cad_cell* cell) {
  if (cell->value) {
    lp_value_destruct(cell->value);
  }

  for (size_t i = 0; i < cell->num_children; i++) {
    cad_tree_destruct(&(cell->children[i]));
  }

  if (cell->children) {
    free(cell->children);
  }
}
