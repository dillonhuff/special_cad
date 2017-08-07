#include "utils.h"

#include <stdlib.h>

lp_polynomial_t* pl_new(const lp_polynomial_context_t* ctx) {
  return lp_polynomial_new(ctx);
}

lp_integer_t mk_int(size_t i) {
  lp_integer_t it;
  lp_integer_construct_from_int(lp_Z, &it, i);

  return it;
}

void print_poly(const lp_polynomial_t* p) {
  lp_polynomial_print(p, stdout);
}


void free_list(void** list, const size_t len) {
  for (size_t i = 0; i < len; i++) {
    free(list[i]);
  }

  free(list);
}

lp_polynomial_t** poly_ptr_list(const size_t len) {
  return (lp_polynomial_t**)(malloc(sizeof(lp_polynomial_t*)*len));
}

void print_poly_list(lp_polynomial_t * const * const ps,
		     const size_t len) {
  for (size_t i = 0; i < len; i++) {
    lp_polynomial_print(ps[i], stdout);
    printf("\n");
  }
}


void pl_delete(pl p) {
  lp_polynomial_delete(p);
}

pl pl_simple_new(const lp_polynomial_context_t* ctx,
		 const size_t a,
		 const lp_variable_t x,
		 const long power) {
  pl p = pl_new(ctx);

  lpint it = mk_int(a);

  lp_polynomial_construct_simple(p, ctx, &it, x, power);

  lp_integer_destruct(&it);
  
  return p;
}


void set_integer_value(lp_assignment_t* asg, const lp_variable_t x, const int i) {
  lpint d_i = mk_int(i);
  lp_value_t* d_value = (lp_value_t*) malloc(sizeof(lp_value_t));
  lp_value_construct_none(d_value);
  lp_value_construct(d_value, LP_VALUE_INTEGER, &d_i);

  lp_assignment_set_value(asg, x, d_value);

  /* lp_integer_destruct(&d_i); */
  /* lp_value_destruct(&d_value); */
  
}
