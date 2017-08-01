#include <stdio.h>

#include <poly/polynomial.h>

int main() {
  lp_variable_db_t var_db;
  /* lp_variable_order_t var_order; */

  lp_polynomial_context_t ctx; //* ctx =
  //lp_polynomial_context_new(lp_Z, &var_db, var_order);
  
  printf("Printing\n");

  lp_polynomial_context_detach(&ctx);
}
