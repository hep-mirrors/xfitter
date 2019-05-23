namespace xfitter{
extern bool isQCDNUMinitialized;
/*!
\brief Initialize QCDNUM if necessary
\details If QCDNUM is not initialized, initializes it by calling QCINIT. Otherwise does nothing.
  Whether QCDNUM has been initalized or not is recorded in global variable "bool isQCDNUMinitialized"
*/
void initQCDNUM();
}
