#include <stdio.h>

#include <poly/polynomial.h>
#include <poly/variable_order.h>
#include <poly/variable_db.h>
#include <poly/polynomial_context.h>
#include <poly/poly.h>

int main() {
  lp_variable_db_t* var_db = lp_variable_db_new();
  lp_variable_order_t* var_order = lp_variable_order_new();

  lp_polynomial_context_t* ctx =
    lp_polynomial_context_new(lp_Z, var_db, var_order);
  
  printf("Printing a number\n");

  lp_integer_t it;
  lp_integer_construct_from_int(lp_Z, &it, 23);

  lp_integer_print(&it, stdout);
  printf("\n");

  lp_polynomial_t* poly = lp_polynomial_new(ctx);
  lp_polynomial_construct_simple(poly, ctx, &it, 1, 1);

  lp_polynomial_print(poly, stdout);
  printf("\n");

  lp_polynomial_context_detach(ctx);
  lp_variable_db_detach(var_db);
  lp_variable_order_detach(var_order);
}
