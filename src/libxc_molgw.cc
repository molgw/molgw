



#include <cstdlib>
#include <cstdio>
#include <xc.h>


extern "C" {

xc_func_type * xc_func_type_malloc()
 { return (xc_func_type *) malloc(sizeof(xc_func_type)); }

void xc_func_type_free(xc_func_type **xc_func)
 { free(*xc_func); *xc_func = NULL; }


}


extern "C" {
int get_family_id(xc_func_type *func)
  { return func->info->family; }
}

extern "C" {
int get_nspin(xc_func_type *func)
  { printf(" from C with love %d spin FBFB\n", func->nspin);
    return func->nspin; }
}

extern "C" {
int set_nspin(xc_func_type *func, int nspin)
  { func->nspin = nspin;
    printf(" from set with love %d %d spin FBFB\n", nspin, func->nspin);
    return 123; }
}
