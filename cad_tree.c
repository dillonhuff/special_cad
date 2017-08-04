#include "cad_tree.h"

#include <assert.h>
#include <stdlib.h>

#include <poly/polynomial.h>

cad_cell make_cad_cell(cad_cell* parent,
		       const size_t num_children,
		       lp_value_t* value) {
  cad_cell c;
  c.parent = parent;
  c.children = (cad_cell*)(malloc(sizeof(cad_cell)*num_children));
  c.num_children = num_children;
  c.value = value;

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
    
    lp_dyadic_rational_t tmp;
    lp_dyadic_rational_construct_from_int(&tmp, 1, 0);

    lp_algebraic_number_t res;
    lp_algebraic_number_construct_from_dyadic_rational(&res, &tmp);

    lp_value_construct(&(test_points[0]), LP_VALUE_ALGEBRAIC, &res);

    lp_dyadic_rational_destruct(&tmp);

    return test_points;
  }

  assert(num_roots > 0);

  *num_test_points_ptr = 2*num_roots + 1;

  lp_value_t* test_points =
    (lp_value_t*)(malloc(sizeof(lp_value_t)*(*num_test_points_ptr)));

  lp_algebraic_number_t neg_inf_pt;
  minus_one(&neg_inf_pt, &(all_roots[0].value.a));

  printf("Minus inf point = ");
  lp_algebraic_number_print(&neg_inf_pt, stdout);
  printf("\n");

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
    
    lp_value_construct(&(test_points[index]), LP_VALUE_ALGEBRAIC, mid_a);
    index++;
  }

  test_points[*num_test_points_ptr - 2] = all_roots[num_roots - 1];

  lp_algebraic_number_t pos_inf_pt;
  plus_one(&pos_inf_pt, &(all_roots[num_roots - 1].value.a));

  printf("Pos inf point = ");
  lp_algebraic_number_print(&pos_inf_pt, stdout);
  printf("\n");

  lp_value_construct(&(test_points[*num_test_points_ptr - 1]), LP_VALUE_ALGEBRAIC, &pos_inf_pt);

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
    return;
  }

  projection_set p = projection_sets[0];

  assert(p.length > 0);
  lp_variable_t next_var = lp_polynomial_top_variable(p.polynomials[0]);

  projection_set* sets_left = projection_sets + 1;

  size_t num_roots = 0;
  lp_value_t* all_roots =
    all_sorted_roots(&num_roots, asg, p.polynomials, p.length);

  printf("# of roots = %zu\n", num_roots);

  size_t num_test_points = 0;
  lp_value_t* all_test_points =
    test_points(&num_test_points, all_roots, num_roots);

  printf("# of test points = %zu\n", num_test_points);

  cad_cell* children =
    (cad_cell*)(malloc(sizeof(cad_cell) * num_test_points));

  for (size_t i = 0; i < num_test_points; i++) {

    assert(all_test_points[i].type == LP_VALUE_ALGEBRAIC);
    lp_algebraic_number_print(&(all_test_points[i].value.a), stdout);
    printf("\n");


    // TODO: Actuall construct new cells
    children[i] = make_cad_cell(root, 0, &(all_test_points[i]));
    lp_assignment_set_value(asg, next_var, &(all_test_points[i]));

    lift_polynomials(&(children[i]), sets_left, num_projection_sets - 1, asg);

    lp_assignment_set_value(asg, next_var, NULL);

    //lp_value_destruct(&fresh_value);
  }

  root->children = children;
  root->num_children = num_test_points;
  
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
    assert(is_algebraic(*(root->value)));

    tab_out(level);
    printf("CELL VALUE = ");
    lp_algebraic_number_print(&((root->value)->value.a), stdout);
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

