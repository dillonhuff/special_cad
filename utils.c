#include "utils.h"

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

