// Dummy file required so that useDynLib(stbench, .registration=TRUE) doesn't fail on empty 'src'

#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>

void attribute_visible R_init_stbench(DllInfo *dll) {
  R_registerRoutines(dll, NULL, NULL, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
