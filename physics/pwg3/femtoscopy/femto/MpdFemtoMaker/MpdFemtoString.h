/**
 * \class MpdFemtoString
 * \brief A file that returns a string
 *
 * The MpdFemtoString header provides a typedef for the string
 *
 * \author Grigory Nigmatkulov (NRNU MEPhI)
 * \date May 18, 2019
 * \email nigmatkulov@gmail.com
 */

#ifndef MpdFemtoString_h
#define MpdFemtoString_h

#ifndef __CINT__

#ifndef MpdFemtoString_noCint
#define MpdFemtoString_noCint

// C++ headers
#include <string>
#if !defined(ST_NO_NAMESPACES)
using std::string;
#endif
typedef string MpdFemtoString;

#endif //MpdFemtoString_noCint

#else // else to #ifndef __CINT__

#ifndef MpdFemtoString_yesCint
#define MpdFemtoString_yesCint
class MpdFemtoString;
#endif //#define MpdFemtoString_yesCint

#endif // #ifndef __CINT__

#endif // #define MpdFemtoString_h
