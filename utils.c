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

