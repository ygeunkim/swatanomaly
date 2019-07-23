#include <Rcpp.h>
using namespace Rcpp;

class som_rect {
private:
  int xgrid;
  int ygrid;
  int num_vars;

public:
  som_rect(int xgrid, int ygrid, int num_vars)
    : xgrid(xgrid), ygrid(ygrid), num_vars(num_vars) {
    // coordinates?
  }

  som_rect(som_rect &copy)
    : xgrid(copy.xgrid), ygrid(copy.ygrid), num_vars(copy.ygrid) {
    // copy init
  }


};
