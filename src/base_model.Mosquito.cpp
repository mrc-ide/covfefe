
#include "base_model.Mosquito.h"

using namespace std;

//------------------------------------------------
// default constructor for mosquito class
Mosquito::Mosquito() {
  host_ID = -1;
  host_infections = -1;
  infection_time = -1;
}

//------------------------------------------------
// constructor for mosquito class
Mosquito::Mosquito(int host_ID, int host_infections, int infection_time) :
  host_ID(host_ID),
  host_infections(host_infections),
  infection_time(infection_time)
{}
