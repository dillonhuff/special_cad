#include "polynomial_examples.h"

#include <poly/variable_db.h>

// A_1, B_1, C_1, D_1, E_1, F_1
lp_polynomial_t** build_2_conic_sections(const lp_polynomial_context_t* ctx,
					 lp_variable_db_t* var_db,
					 lp_variable_order_t* var_order) {

  // Create variables
  lp_variable_t x = lp_variable_db_new_variable(var_db, "x");
  lp_variable_t y = lp_variable_db_new_variable(var_db, "y");

  lp_variable_t A_1 = lp_variable_db_new_variable(var_db, "A_1");

  lp_variable_order_push(var_order, A_1);
  lp_variable_order_push(var_order, x);  
  lp_variable_order_push(var_order, y);

  
}
