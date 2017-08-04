#pragma once

#include <poly/value.h>

#include "utils.h"

struct cad_cell {
  struct cad_cell* parent;
  struct cad_cell* children;
  size_t num_children;

  lp_value_t* value;
};

typedef struct cad_cell cad_cell;

struct projection_set {
  size_t length;
  lp_polynomial_t** polynomials;
};

typedef struct projection_set projection_set;

cad_cell make_cad_cell(cad_cell* parent,
		       const size_t num_children,
		       lp_value_t* value);

projection_set make_projection_set(lp_polynomial_t** polys,
				   size_t length);

void lift_polynomials(cad_cell* root,
		      projection_set* projection_sets,
		      const size_t num_projection_sets,
		      lp_assignment_t* asg);

void print_cad_tree(cad_cell* root);

size_t is_algebraic(const lp_value_t val);

lp_value_t* all_sorted_roots(size_t* num_roots_ptr,
			     lp_assignment_t* const assignment,
			     lp_polynomial_t * const * const ps,
			     const size_t ps_len);

lp_value_t* test_points(size_t* num_test_points_ptr,
			lp_value_t const * const all_roots,
			const size_t num_roots);
