// -------------------------------------------------------------------------
// -----          TShieldGeometry source file                          -----
// -----          Created by D. Sosnov                                 -----
// -------------------------------------------------------------------------

#include "TShieldGeometry.h"
namespace tgeanttoshield {
// Operators for comparsion floating point values
bool doubleNE(const double &l, const double &r)
{
   return (abs(l - r) > 1E-15);
}
bool doubleEQ(const double &l, const double &r)
{
   return !doubleNE(l, r);
}
bool doubleLT(const double &l, const double &r)
{
   return (doubleNE(l, r) && (l < r));
}
bool doubleLE(const double &l, const double &r)
{
   return ((l < r) || doubleEQ(l, r));
}
bool doubleGT(const double &l, const double &r)
{
   return !doubleLE(l, r);
}
bool doubleGE(const double &l, const double &r)
{
   return !doubleLT(l, r);
}

TGeoHMatrix operator*(TGeoHMatrix &matrixR, const TGeoMatrix &matrixL)
{
   TGeoHMatrix mm = TGeoHMatrix(matrixR);
   mm.Multiply(&matrixL);
   return TGeoHMatrix(mm);
}
TGeoHMatrix operator*(TGeoHMatrix &matrixR, const TGeoHMatrix &matrixL)
{
   TGeoHMatrix mm = TGeoHMatrix(matrixR);
   return mm *= matrixL;
}
TGeoHMatrix operator*(const TGeoMatrix &matrixR, const TGeoMatrix &matrixL)
{
   TGeoHMatrix hmr = TGeoHMatrix(matrixR);
   return hmr * matrixL;
}
TGeoHMatrix operator*(const TGeoMatrix &matrixR, const TGeoHMatrix &matrixL)
{
   TGeoHMatrix hmr = TGeoHMatrix(matrixR);
   return hmr * matrixL;
}
TGeoTranslation operator+(const TGeoTranslation &matrixR, const TGeoTranslation &matrixL)
{
   TGeoTranslation m = TGeoTranslation(matrixR);
   m.Add(&matrixL);
   return TGeoTranslation(m);
}
TGeoTranslation operator-(const TGeoTranslation matrixR, const TGeoTranslation &matrixL)
{
   return matrixR + matrixL.Inverse();
}
TGeoTranslation operator+(const TGeoTranslation matrixR, const TGeoMatrix &matrixL)
{
   return matrixR + TGeoTranslation(matrixL);
}
TGeoTranslation operator*(const TGeoTranslation matrixR, const Double_t &scale)
{
   Double_t *tmp1 = (Double_t *)matrixR.GetTranslation();
   return TGeoTranslation(tmp1[0] * scale, tmp1[1] * scale, tmp1[2] * scale);
}
void addVectorToElement(SGeoBody &body, unsigned int position, TGeoTranslation vector)
{
   body.parameters[position]     = (double)(vector.GetTranslation()[0]);
   body.parameters[position + 1] = (double)(vector.GetTranslation()[1]);
   body.parameters[position + 2] = (double)(vector.GetTranslation()[2]);
}
} // namespace tgeanttoshield