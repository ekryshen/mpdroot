#ifndef __INTERVAL_TREE_H
#define __INTERVAL_TREE_H

//------------------------------------------------------------------------------------------------------------------------
/// \class IntervalTree
///
/// \brief
/// \author Sergei Lobastov (LHE, JINR, Dubna)
//------------------------------------------------------------------------------------------------------------------------

#include <assert.h>
#include <vector>
#include <algorithm>
#include <iostream>
//------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------
//		class	Interval
//------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------
template <class T, typename K = double>
class Interval {
public:
   K start;
   K stop;
   T value;

   Interval(K s, K e, const T &v) : start(s), stop(e), value(v) {}
};
//------------------------------------------------------------------------------------------------------------------------
template <class T, typename K>
K intervalStart(const Interval<T, K> &i)
{
   return i.start;
}
//------------------------------------------------------------------------------------------------------------------------
template <class T, typename K>
K intervalStop(const Interval<T, K> &i)
{
   return i.stop;
}
//------------------------------------------------------------------------------------------------------------------------
template <class T, typename K>
std::ostream &operator<<(std::ostream &out, const Interval<T, K> &i)
{
   out << "Interval(" << i.start << ", " << i.stop << "): " << i.value;
   return out;
}
//------------------------------------------------------------------------------------------------------------------------
template <class T, typename K = double>
class IntervalStartSorter {
public:
   bool operator()(const Interval<T, K> &a, const Interval<T, K> &b) { return a.start < b.start; }
};
//------------------------------------------------------------------------------------------------------------------------
template <class T, typename K = double>
class IntervalValueSorter {
public:
   bool operator()(const Interval<T, K> &a, const Interval<T, K> &b) { return a.value < b.value; }
};
//------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------
//		class	IntervalTree
//------------------------------------------------------------------------------------------------------------------------
//------------------------------------------------------------------------------------------------------------------------
template <class T, typename K = double>
class IntervalTree {
public:
   typedef Interval<T, K>        interval;
   typedef std::vector<interval> intervalVector;
   typedef IntervalTree<T, K>    intervalTree;

   intervalVector intervals;
   intervalTree  *left;
   intervalTree  *right;
   K              center;

   //------------------------------------------------------------------------------------------------------------------------
   IntervalTree<T, K>(void) : left(nullptr), right(nullptr), center(0.) {}
   //------------------------------------------------------------------------------------------------------------------------
   IntervalTree<T, K>(const intervalTree &other) : left(nullptr), right(nullptr)
   {
      center    = other.center;
      intervals = other.intervals;
      if (other.left) left = new intervalTree(*other.left);
      if (other.right) right = new intervalTree(*other.right);
   }
   //------------------------------------------------------------------------------------------------------------------------
   ~IntervalTree(void)
   {
      delete left;
      delete right;
   }
   //------------------------------------------------------------------------------------------------------------------------
   IntervalTree<T, K> &operator=(const intervalTree &other)
   {
      center    = other.center;
      intervals = other.intervals;

      if (other.left)
         left = new intervalTree(*other.left);
      else {
         if (left) delete left;
         left = nullptr;
      }

      if (other.right)
         right = new intervalTree(*other.right);
      else {
         if (right) delete right;
         right = nullptr;
      }

      return *this;
   }
   //------------------------------------------------------------------------------------------------------------------------
   IntervalTree<T, K>(intervalVector &ivals, std::size_t depth = 16, std::size_t minbucket = 64, K leftextent = 0,
                      K rightextent = 0, std::size_t maxbucket = 512)
      : left(nullptr), right(nullptr)
   {
      --depth;
      IntervalStartSorter<T, K> intervalStartSorter;
      if (depth == 0 || (ivals.size() < minbucket && ivals.size() < maxbucket)) {
         std::sort(ivals.begin(), ivals.end(), intervalStartSorter);
         intervals = ivals;
      } else {
         if (leftextent == 0 && rightextent == 0) {
            // sort intervals by start
            std::sort(ivals.begin(), ivals.end(), intervalStartSorter);
         }

         K leftp   = 0.;
         K rightp  = 0.;
         K centerp = 0.;

         if (leftextent || rightextent) {
            leftp  = leftextent;
            rightp = rightextent;
         } else {
            leftp = ivals.front().start;
            std::vector<K> stops;
            stops.resize(ivals.size());
            transform(ivals.begin(), ivals.end(), stops.begin(), intervalStop<T, K>);
            rightp = *max_element(stops.begin(), stops.end());
         }

         // centerp = ( leftp + rightp ) / 2;
         centerp = ivals.at(ivals.size() / 2).start;
         center  = centerp;

         intervalVector lefts;
         intervalVector rights;

         for (typename intervalVector::iterator it = ivals.begin(); it != ivals.end(); ++it) {
            interval &range = *it;

            if (range.stop < center)
               lefts.push_back(range);
            else if (range.start > center)
               rights.push_back(range);
            else
               intervals.push_back(range);
         }

         if (!lefts.empty()) left = new intervalTree(lefts, depth, minbucket, leftp, centerp);
         if (!rights.empty()) right = new intervalTree(rights, depth, minbucket, centerp, rightp);
      }
   }
   //------------------------------------------------------------------------------------------------------------------------
   void findOverlapping(K start, K stop, intervalVector &overlapping) const
   {
      if (!intervals.empty() && !(stop < intervals.front().start)) {
         for (typename intervalVector::const_iterator it = intervals.begin(); it != intervals.end(); ++it) {
            const interval &range = *it;
            if (range.stop >= start && range.start <= stop) overlapping.push_back(range);
         }
      }

      if (left && start <= center) left->findOverlapping(start, stop, overlapping);
      if (right && stop >= center) right->findOverlapping(start, stop, overlapping);
   }
   //------------------------------------------------------------------------------------------------------------------------
   void findContained(K start, K stop, intervalVector &contained) const
   {
      if (!intervals.empty() && !(stop < intervals.front().start)) {
         for (typename intervalVector::const_iterator it = intervals.begin(); it != intervals.end(); ++it) {
            const interval &range = *it;
            if (range.start >= start && range.stop <= stop) contained.push_back(range);
         }
      }

      if (left && start <= center) left->findContained(start, stop, contained);
      if (right && stop >= center) right->findContained(start, stop, contained);
   }
   //------------------------------------------------------------------------------------------------------------------------
   void dump(const char *comment = nullptr, std::ostream &out = std::cout, int treeLevel = 0) const
   {
      treeLevel++;

      if (comment) out << comment;
      out << " left= " << left << " right= " << right << "  size= " << intervals.size();
      for (typename intervalVector::const_iterator it = intervals.begin(), itEnd = intervals.end(); it != itEnd; it++)
         out << std::endl << (*it);

      if (left != nullptr) {
         out << "\n Left_" << treeLevel;
         left->dump(comment, out);
      }

      if (right != nullptr) {
         out << "\n Right_" << treeLevel;
         right->dump(comment, out);
      }
   }
   //------------------------------------------------------------------------------------------------------------------------
   static void dumpIntervals(const intervalVector *pVector, const char *comment = nullptr,
                             std::ostream &out = std::cout)
   {
      assert(pVector != nullptr);

      if (comment) out << comment;
      out << "  size= " << pVector->size();
      for (typename intervalVector::const_iterator it = pVector->begin(), itEnd = pVector->end(); it != itEnd; it++)
         out << std::endl << (*it);
   }
   //------------------------------------------------------------------------------------------------------------------------
};
//------------------------------------------------------------------------------------------------------------------------
#endif
