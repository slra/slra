#include <limits>
#include <memory.h>
#include <cstdarg>
extern "C" {
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>
}
#include "slra.h"

typedef Structure* pStructure;

pStructure * MosaicHStructure::allocStripe( gsl_vector *m_l, gsl_vector *n_k, 
                 gsl_vector *w )  {
  pStructure *res = new pStructure[n_k->size];

  for (size_t k = 0; k < n_k->size; k++) {
    res[k] = new  LayeredHStructure(m_l->data, m_l->size, n_k->data[k], 
                                   (w != NULL ? w->data : NULL));
    if (w != NULL && (w->size != m_l->size)) {
      w->data += m_l->size;
    }
  }
  return res;
}
 
MosaicHStructure:: MosaicHStructure( gsl_vector *m_l, gsl_vector *n_k,  
     gsl_vector *w ) : StripedStructure(n_k->size, allocStripe(m_l, n_k, w),
         w == NULL || w->size == m_l->size)  {
}


