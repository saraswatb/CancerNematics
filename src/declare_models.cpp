// ======================================================================
// Model declaration (gets compiled in models.cpp)
// model headers and declare_model must be consistent!


// model headers
#include "models/ViscLyo.hpp"

void DeclareModels()
{
  declare_model<ViscLyo>(
      "ViscLyo",
      "Lyotropic model (see Biphasic, lyotropic, active nematics. PRL Blow et al. 2015) - but the passive phase is viscoelastic."      
      );
}
