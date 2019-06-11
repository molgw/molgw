



#include <cstdlib>
#include <xc.h>


extern "C" {

xc_func_type * xc_func_type_malloc()
 { return (xc_func_type *) malloc(sizeof(xc_func_type)); }
 
void xc_func_type_free(xc_func_type **xc_func)
 { free(*xc_func); *xc_func = NULL; }


}


