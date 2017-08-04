#include "cad_tree.h"

#include <assert.h>
#include <stdlib.h>

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

void lift_polynomials(cad_cell* root,
		      projection_set* projection_sets,
		      const size_t num_projection_sets,
		      lp_assignment_t* asg) {
  
}

void print_cad_tree(cad_cell* root) {
  if (root->value == NULL) {
    printf("ROOT\n");
  } else {
    assert(is_algebraic(*(root->value)));
    lp_algebraic_number_print(&((root->value)->value.a), stdout);
    printf("\n");
  }

  for (size_t i = 0; i < root->num_children; i++) {
    print_cad_tree(&(root->children[i]));
  }
}

size_t is_algebraic(const lp_value_t val) {
  return val.type == LP_VALUE_ALGEBRAIC;
}

