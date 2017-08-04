#include "cad_tree.h"

#include <stdlib.h>

cad_cell make_cad_cell(cad_cell* parent,
		       const size_t num_children,
		       lp_value_t* value) {
  cad_cell c;
  c.parent = parent;
  c.children = (cad_cell*)(malloc(sizeof(cad_cell)*num_children));
  c.value = value;

  return c;
}
