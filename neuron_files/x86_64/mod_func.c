#include <stdio.h>
#include "hocdec.h"
extern int nrnmpi_myid;
extern int nrn_nobanner_;
modl_reg(){
  if (!nrn_nobanner_) if (nrnmpi_myid < 1) {
    fprintf(stderr, "Additional mechanisms from files\n");

    fprintf(stderr," AType_potassium.mod");
    fprintf(stderr," HH_traub.mod");
    fprintf(stderr," IL_gutnick.mod");
    fprintf(stderr," IM_cortex.mod");
    fprintf(stderr," cadecay.mod");
    fprintf(stderr," cadecay_destexhe.mod");
    fprintf(stderr," cagk2.mod");
    fprintf(stderr," flushf.mod");
    fprintf(stderr," intfire4r.mod");
    fprintf(stderr," kA.mod");
    fprintf(stderr," kM.mod");
    fprintf(stderr," kca.mod");
    fprintf(stderr," kfasttab.mod");
    fprintf(stderr," kslowtab.mod");
    fprintf(stderr," lcafixed.mod");
    fprintf(stderr," nafast.mod");
    fprintf(stderr," nagran.mod");
    fprintf(stderr," nmdanet.mod");
    fprintf(stderr," odorinput.mod");
    fprintf(stderr, "\n");
  }
  _AType_potassium_reg();
  _HH_traub_reg();
  _IL_gutnick_reg();
  _IM_cortex_reg();
  _cadecay_reg();
  _cadecay_destexhe_reg();
  _cagk2_reg();
  _flushf_reg();
  _intfire4r_reg();
  _kA_reg();
  _kM_reg();
  _kca_reg();
  _kfasttab_reg();
  _kslowtab_reg();
  _lcafixed_reg();
  _nafast_reg();
  _nagran_reg();
  _nmdanet_reg();
  _odorinput_reg();
}
