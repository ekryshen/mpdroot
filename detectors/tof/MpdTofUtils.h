//--------------------------------------------------------------------------------------------------------------
#ifndef __HH_MPDTOFUTILS
#define __HH_MPDTOFUTILS

//#include <limits>
#include <Rtypes.h> // bit manipulation
//--------------------------------------------------------------------------------------------------------------
namespace MpdTofUtils {
enum k_side { Up = 0, Right, Down, Left, Absent = -1 };
enum k_HitType {
   IsSingle   = BIT(0),
   IsDouble   = BIT(1),
   HaveTail   = BIT(3),
   InCluster  = BIT(4),
   McAbsent   = BIT(5),
   IsSelected = BIT(12)
};
enum k_LinkType { mcPointIndex = 1, mcTrackIndex = 2, volumeUID = 3, expDigitIndex = 4 };
enum k_DetType { IsBarrel = 1, IsEndcapLeft = 2, IsEndcapRight = 3 };

//	constexpr Double_t nan = std::numeric_limits<double>::quiet_NaN();
};     // namespace MpdTofUtils
#endif // #ifndef __HH_MPDTOFUTILS
