#include <limits>
#include <memory.h>
#include <math.h>
extern "C" {
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_math.h>
}
#include "slra.h"

typedef Structure* pStructure;

pStructure *WMosaicHStructure::allocStripe( gsl_vector *m_l, gsl_vector *n_k,  
                 gsl_vector *w )  {
  pStructure *res = new pStructure[n_k->size];

  for (size_t k = 0; k < n_k->size; k++) {
    res[k] = new WLayeredHStructure(m_l->data, m_l->size, n_k->data[k], 
                                    w);
    if (w != NULL) {
      w->data += res[k]->getNp();
    }
  }
  return res;
}

WMosaicHStructure::WMosaicHStructure( gsl_vector *m_l, gsl_vector *n_k,  
    gsl_vector *w ) : StripedStructure(n_k->size, allocStripe(m_l, n_k, w)) {
    
}
