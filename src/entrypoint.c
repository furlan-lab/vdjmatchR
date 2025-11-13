// Forward routine registration from C to the extendr-generated function.
void R_init_vdjmatchR_extendr(void *dll);

void R_init_vdjmatchR(void *dll) {
    R_init_vdjmatchR_extendr(dll);
}
