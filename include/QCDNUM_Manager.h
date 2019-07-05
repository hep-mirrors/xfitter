namespace xfitter{
extern bool isQCDNUMinitialized;
/*!
\brief Initialize QCDNUM if necessary
\details If QCDNUM is not initialized, initializes it by calling QCINIT. Otherwise does nothing.
  Whether QCDNUM has been initalized or not is recorded in global variable "bool isQCDNUMinitialized"
*/
void initQCDNUM();

extern bool isZMSTFinitialized;
/*!
\brief Initialize ZMSTF if necessary
\details ZMSTF is a QCDNUM addon for calculation of structure functions.
  If ZMSTF is not initialized, initializes it by calling ZMFILLW. Otherwise does nothing.
  Calls initQCDNUM if necessary.
  Whether ZMSTF has been initalized or not is recorded in global variable "bool isZMSTFinitialized"
*/
void requireZMSTF();

}
