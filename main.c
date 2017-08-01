#include <stdio.h>

#include <poly/variable_order.h>
#include <poly/variable_db.h>
#include <poly/polynomial_context.h>
#include <poly/poly.h>



int main() {
  lp_variable_db_t* var_db = lp_variable_db_new();
  lp_variable_order_t* var_order = lp_variable_order_new();

  lp_polynomial_context_t* ctx =
    lp_polynomial_context_new(lp_Z, var_db, var_order);
  
  printf("Printing\n");

  lp_polynomial_context_detach(ctx);
  lp_variable_db_detach(var_db);
}
