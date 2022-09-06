#include <cstdarg>
#include <cstdio>

#include "Debug.h"

namespace TpcAlignment {

Debug::Debug(bool mode) : vDebugMode(mode) {}

Debug::~Debug() {}

void Debug::Print(char const *const _Format, ...) const
{
   if (!vDebugMode) {
      return;
   }
   va_list args;
   va_start(args, _Format);
   vprintf(_Format, args);
   va_end(args);
}
} // namespace TpcAlignment
