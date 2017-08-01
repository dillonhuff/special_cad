#include <stdio.h>
#include <stdlib.h>

#include <poly/assignment.h>
#include <poly/polynomial.h>
#include <poly/variable_order.h>
#include <poly/variable_db.h>
#include <poly/polynomial_context.h>
#include <poly/upolynomial.h>
#include <poly/poly.h>

int main() {
  lp_variable_db_t* var_db = lp_variable_db_new();

  // Create variables
  lp_variable_t u = lp_variable_db_new_variable(var_db, "u");  

  // Create variable order
  lp_variable_order_t* var_order = lp_variable_order_new();
  lp_variable_order_push(var_order, u);

  lp_polynomial_context_t* ctx =
    lp_polynomial_context_new(lp_Z, var_db, var_order);
  
  printf("Printing a number\n");

  lp_integer_t it;
  lp_integer_construct_from_int(lp_Z, &it, 23);

  lp_integer_t one;
  lp_integer_construct_from_int(lp_Z, &one, 1);
  
  lp_integer_print(&it, stdout);
  printf("\n");

  // Build up polynomial 23u^2 + u from monomials
  /* lp_polynomial_t* x2 = lp_polynomial_new(ctx); */
  /* lp_polynomial_construct_simple(x2, ctx, &it, u, 2); */

  /* lp_polynomial_print(x2, stdout); */
  /* printf("\n"); */
  
  /* lp_polynomial_t* x = lp_polynomial_new(ctx); */
  /* lp_polynomial_construct_simple(x, ctx, &one, u, 1); */

  /* lp_polynomial_print(x, stdout); */
  /* printf("\n"); */
  
  /* lp_polynomial_t* poly = lp_polynomial_new(ctx); */
  /* lp_polynomial_add(poly, x2, x); */

  /* lp_polynomial_print(poly, stdout); */
  /* printf("\n"); */

  lp_upolynomial_t* x2 = lp_upolynomial_construct_power(lp_Z, 23, 2);
  lp_upolynomial_t* x = lp_upolynomial_construct_power(lp_Z, 1, 1);

  lp_upolynomial_t* poly = lp_upolynomial_add(x2, x);
  size_t deg = lp_upolynomial_degree(poly);

  printf("Degree of poly = %zu\n", deg);

  // Isolate roots
  lp_assignment_t* assignment = lp_assignment_new(var_db);

  /* void* data = malloc(sizeof(lp_rational_t)*deg); */
  /* lp_value_t* roots = lp_value_new(LP_VALUE_RATIONAL, data); */

  lp_algebraic_number_t* roots =
    (lp_algebraic_number_t*)(malloc(sizeof(lp_algebraic_number_t)*deg));

  size_t roots_size;
  //lp_polynomial_roots_isolate(poly, assignment, roots, &roots_size);
  lp_upolynomial_roots_isolate(poly, roots, &roots_size);

  printf("# of roots of poly = %zu\n", roots_size);

  for (size_t i = 0; i < roots_size; i++) {
    printf("ROOT\n");
    lp_algebraic_number_print(&(roots[i]), stdout);
    //printf("root type = %d\n", roots[i].type);
    //printf("root approximation = %f\n", lp_value_to_double(&(roots[i])));
  }

  // Cleanup
  
  lp_assignment_delete(assignment);
  lp_polynomial_context_detach(ctx);
  lp_variable_db_detach(var_db);
  lp_variable_order_detach(var_order);


}
