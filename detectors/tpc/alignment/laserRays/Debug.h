#ifndef DEBUG_HH
#define DEBUG_HH

namespace TpcAlignmentLaserRays {
/// <summary>
/// Output Debug info
/// </summary>
class Debug {
private:
   bool vDebugMode{false};

public:
   /// <summary>
   /// If mode is true then ouput debug data.
   /// </summary>
   /// <param name="mode">Switch on/off debug output</param>
   Debug(bool mode);
   ~Debug();
   /// <summary>
   /// If true then print all messages
   /// </summary>
   /// <param name="mode"></param>
   void SetDebugMode(bool mode);

   /// <summary>
   /// Prints debug message if allowed.
   /// </summary>
   /// <param name="_Format">same as for printf function</param>
   /// <param name="">Arguments for printf function</param>
   void Print(char const *const _Format, ...) const;
};
} // namespace TpcAlignmentLaserRays
#endif
