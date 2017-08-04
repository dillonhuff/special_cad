#pragma once

struct cad_cell {
  struct cad_cell* parent;
  struct cad_cell* children;

  
};

typedef struct cad_cell cad_cell;
