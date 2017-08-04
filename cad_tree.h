#pragma once

#include <poly/value.h>

struct cad_cell {
  struct cad_cell* parent;
  struct cad_cell* children;

  lp_value_t value;
};

typedef struct cad_cell cad_cell;
