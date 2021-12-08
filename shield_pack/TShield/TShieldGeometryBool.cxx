// -------------------------------------------------------------------------
// -----          TShieldGeometry source file                          -----
// -----          Created by D. Sosnov                                 -----
// -------------------------------------------------------------------------

#include "TShieldGeometry.h"
namespace tgeanttoshield {
zoneElement notElement(zoneElement el)
{
   return std::pair<int, SGeoBody>(el.first * -1, el.second);
}
zoneList notZone(zoneList list)
{
   int cur, count = 1;
   for (unsigned int i = 0; i < list.size(); ++i) {
      count *= list.at(i).size();
   }
   zoneList                 out;
   std::vector<zoneElement> tmp;
   for (int k = 0; k < count; ++k) {
      cur = k;
      tmp.clear();
      for (unsigned int i = 0; i < list.size(); ++i) {
         int cursize = list.at(i).size();
         tmp.push_back(notElement(list.at(i).at(cur % cursize)));
         cur /= cursize;
      }
      out.push_back(tmp);
   }
   return out;
}
zoneList orZone(zoneList list1, zoneList list2)
{
   zoneList outList = list1;
   outList.insert(outList.end(), list2.begin(), list2.end());
   return outList;
}
zoneList andZone(zoneList list1, zoneList list2)
{
   zoneList                 outList;
   std::vector<zoneElement> tmp;
   for (unsigned int i = 0; i < list1.size(); ++i) {
      for (unsigned int j = 0; j < list2.size(); ++j) {
         tmp = list1.at(i);
         tmp.insert(tmp.end(), list2.at(j).begin(), list2.at(j).end());
         outList.push_back(tmp);
      }
   }
   return outList;
}
} // namespace tgeanttoshield