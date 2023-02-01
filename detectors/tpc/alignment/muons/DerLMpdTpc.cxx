/// \class DerLMpdTpc
/// \brief Partial derivatives for TPC alignment
///
/// \author Valentin Kuzmin, SINP of Moscow State University

#include <iostream>
using namespace std;

#include <TMath.h>
#include "DerLMpdTpc.h"

void DerLMpdTpc::Init()
{

   Double_t c1, s1, c2, s2, c3, s3, tl2g[3][3];

   for (Int_t is = 0; is < 24; is++) {
      TVector3 r0_loc;
      Double_t alpha, beta, gamma;
      fTpcSecGeo->GetAlignment(is, r0_loc, alpha, beta, gamma);

      c1 = TMath::Cos(alpha);
      s1 = TMath::Sin(alpha);
      c2 = TMath::Cos(beta);
      s2 = TMath::Sin(beta);
      c3 = TMath::Cos(gamma);
      s3 = TMath::Sin(gamma);
      fTpcSecGeo->SectorBackTransformation(is, G0_gsc[is], glo2loc[is]);
      fTpcSecGeo->GCS2TLCStransformation(is, G0_tlsc[is], Tglo2lo[is]);
      /*printf("G0_tlsc(%f %f %f)\n",G0_tlsc[is].X(),G0_tlsc[is].Y(),G0_tlsc[is].Z());
      printf("Tglo2lo    %f %f %f\n",Tglo2lo[is](0,0),Tglo2lo[is](0,1),Tglo2lo[is](0,2));
      printf("           %f %f %f\n",Tglo2lo[is](1,0),Tglo2lo[is](1,1),Tglo2lo[is](1,2));
      printf("           %f %f %f\n",Tglo2lo[is](2,0),Tglo2lo[is](2,1),Tglo2lo[is](2,2));*/
      fTpcSecGeo->DesignSectorTransformation(is, Tlo2glo[is]);
      // TRotation dlo2glo is the transformation theoretical sector TLCS -> GCS
      // matrix tl2g=TRotation dlo2glo
      for (Int_t i = 0; i < 3; i++)
         for (Int_t j = 0; j < 3; j++) tl2g[i][j] = Tlo2glo[is](i, j);

      // partial derivatives of the transformation
      //  ||R|| - transformation matrix (LCS->GCS) with alignment
      //  ||tl2g|| - transformation matrix (TLCS->GCS) without alignment
      //  ||dA|| - matrix of derivatives of the alignmet matrix (LCS->TLCS)
      //  ||dR|| - matrix of derivatives of ||dR|| = ||R|| ||dl2g|| (LCS->GCS)

      // r_g=r0_sec+||dl2g||r_is
      // r_is=r0_a+||R^-1||r_s
      // r_g=r0_sec+||dl2g||(r0_a+||R^-1||r_s)
      // r_g=(r0_sec+||dl2g||r0_a)+||dl2g||||R^-1||r_s

      // = = = = alpha dR=dR^(-1)/dp4
      dAdp[0][0][0][is] = -s1 * c2;
      dAdp[0][1][0][is] = -s1 * s2 * s3 - c1 * c3;
      dAdp[0][2][0][is] = -s1 * s2 * c3 + c1 * s3;

      dAdp[1][0][0][is] = c1 * c2;
      dAdp[1][1][0][is] = c1 * s2 * s3 - s1 * c3;
      dAdp[1][2][0][is] = c1 * s2 * c3 + s1 * s3;

      dAdp[2][0][0][is] = 0;
      dAdp[2][1][0][is] = 0;
      dAdp[2][2][0][is] = 0;

      for (Int_t i = 0; i < 3; i++)
         for (Int_t j = 0; j < 3; j++) {
            dRdp[i][j][0][is] = 0;
            for (Int_t k = 0; k < 3; k++) dRdp[i][j][0][is] += tl2g[i][k] * dAdp[k][j][0][is];
         }

      // = = = =  dR=d2||R^(-1)||/dp4dp4
      d2Adp2[0][0][0][0][is] = -c1 * c2;
      d2Adp2[0][1][0][0][is] = -c1 * s2 * s3 + s1 * c3;
      d2Adp2[0][2][0][0][is] = -c1 * s2 * c3 - s1 * s3;

      d2Adp2[1][0][0][0][is] = -s1 * c2;
      d2Adp2[1][1][0][0][is] = -s1 * s2 * s3 - c1 * c3;
      d2Adp2[1][2][0][0][is] = -s1 * s2 * c3 + c1 * s3;

      d2Adp2[2][0][0][0][is] = 0;
      d2Adp2[2][1][0][0][is] = 0;
      d2Adp2[2][2][0][0][is] = 0;

      for (Int_t i = 0; i < 3; i++)
         for (Int_t j = 0; j < 3; j++) {
            d2Rdp2[i][j][0][0][is] = 0;
            for (Int_t k = 0; k < 3; k++) d2Rdp2[i][j][0][0][is] += tl2g[i][k] * d2Adp2[k][j][0][0][is];
         }

      // = = = = dR=d2||R^(-1)||/dp4dp5
      d2Adp2[0][0][0][1][is] = s1 * s2;
      d2Adp2[0][1][0][1][is] = -s1 * c2 * s3;
      d2Adp2[0][2][0][1][is] = -s1 * c2 * c3;

      d2Adp2[1][0][0][1][is] = -c1 * s2;
      d2Adp2[1][1][0][1][is] = c1 * c2 * s3;
      d2Adp2[1][2][0][1][is] = c1 * c2 * c3;

      d2Adp2[2][0][0][1][is] = 0;
      d2Adp2[2][1][0][1][is] = 0;
      d2Adp2[2][2][0][1][is] = 0;

      for (Int_t i = 0; i < 3; i++)
         for (Int_t j = 0; j < 3; j++) {
            d2Rdp2[i][j][0][1][is] = 0;
            for (Int_t k = 0; k < 3; k++) d2Rdp2[i][j][0][1][is] += tl2g[i][k] * d2Adp2[k][j][0][1][is];
            d2Rdp2[i][j][1][0][is] = d2Rdp2[i][j][0][1][is];
            d2Adp2[i][j][1][0][is] = d2Adp2[i][j][0][1][is];
         }

      // = = = = dR=d2||R^(-1)||/dp4dp6
      d2Adp2[0][0][0][2][is] = 0;
      d2Adp2[0][1][0][2][is] = -s1 * s2 * c3 + c1 * s3;
      d2Adp2[0][2][0][2][is] = s1 * s2 * s3 + c1 * c3;

      d2Adp2[1][0][0][2][is] = 0;
      d2Adp2[1][1][0][2][is] = c1 * s2 * c3 + s1 * s3;
      d2Adp2[1][2][0][2][is] = -c1 * s2 * s3 + s1 * c3;

      d2Adp2[2][0][0][2][is] = 0;
      d2Adp2[2][1][0][2][is] = 0;
      d2Adp2[2][2][0][2][is] = 0;

      for (Int_t i = 0; i < 3; i++)
         for (Int_t j = 0; j < 3; j++) {
            d2Rdp2[i][j][0][2][is] = 0;
            for (Int_t k = 0; k < 3; k++) d2Rdp2[i][j][0][2][is] += tl2g[i][k] * d2Adp2[k][j][0][2][is];
            d2Rdp2[i][j][2][0][is] = d2Rdp2[i][j][0][2][is];
            d2Adp2[i][j][2][0][is] = d2Adp2[i][j][0][2][is];
         }

      // = = = = beta dR=dR^(-1)/dp5
      dAdp[0][0][1][is] = -c1 * s2;
      dAdp[0][1][1][is] = c1 * c2 * s3;
      dAdp[0][2][1][is] = c1 * c2 * c3;

      dAdp[1][0][1][is] = -s1 * s2;
      dAdp[1][1][1][is] = s1 * c2 * s3;
      dAdp[1][2][1][is] = s1 * c2 * c3;

      dAdp[2][0][1][is] = -c2;
      dAdp[2][1][1][is] = -s2 * s3;
      dAdp[2][2][1][is] = -s2 * c3;

      for (Int_t i = 0; i < 3; i++)
         for (Int_t j = 0; j < 3; j++) {
            dRdp[i][j][1][is] = 0;
            for (Int_t k = 0; k < 3; k++) dRdp[i][j][1][is] += tl2g[i][k] * dAdp[k][j][1][is];
         }

      // = = = = dR=dR^(-1)/dp5dp5
      d2Adp2[0][0][1][1][is] = -c1 * c2;
      d2Adp2[0][1][1][1][is] = -c1 * s2 * s3;
      d2Adp2[0][2][1][1][is] = -c1 * s2 * c3;

      d2Adp2[1][0][1][1][is] = -s1 * c2;
      d2Adp2[1][1][1][1][is] = -s1 * s2 * s3;
      d2Adp2[1][2][1][1][is] = -s1 * s2 * c3;

      d2Adp2[2][0][1][1][is] = s2;
      d2Adp2[2][1][1][1][is] = -c2 * s3;
      d2Adp2[2][2][1][1][is] = -c2 * c3;

      for (Int_t i = 0; i < 3; i++)
         for (Int_t j = 0; j < 3; j++) {
            d2Rdp2[i][j][1][1][is] = 0;
            for (Int_t k = 0; k < 3; k++) d2Rdp2[i][j][1][1][is] += tl2g[i][k] * d2Adp2[k][j][1][1][is];
         }

      // = = = = dR=dR^(-1)/dp5dp6
      d2Adp2[0][0][1][2][is] = 0;
      d2Adp2[0][1][1][2][is] = c1 * c2 * c3;
      d2Adp2[0][2][1][2][is] = -c1 * c2 * s3;

      d2Adp2[1][0][1][2][is] = 0;
      d2Adp2[1][1][1][2][is] = s1 * c2 * c3;
      d2Adp2[1][2][1][2][is] = -s1 * c2 * s3;

      d2Adp2[2][0][1][2][is] = 0;
      d2Adp2[2][1][1][2][is] = -s2 * c3;
      d2Adp2[2][2][1][2][is] = s2 * s3;

      for (Int_t i = 0; i < 3; i++)
         for (Int_t j = 0; j < 3; j++) {
            d2Rdp2[i][j][1][2][is] = 0;
            for (Int_t k = 0; k < 3; k++) d2Rdp2[i][j][1][2][is] += tl2g[i][k] * d2Adp2[k][j][1][2][is];
            d2Rdp2[i][j][2][1][is] = d2Rdp2[i][j][1][2][is];
            d2Adp2[i][j][2][1][is] = d2Adp2[i][j][1][2][is];
         }

      // = = = = gamma dR=dR^(-1)/dp6
      dAdp[0][0][2][is] = 0;
      dAdp[0][1][2][is] = c1 * s2 * c3 + s1 * s3;
      dAdp[0][2][2][is] = -c1 * s2 * s3 + s1 * c3;

      dAdp[1][0][2][is] = 0;
      dAdp[1][1][2][is] = s1 * s2 * c3 - c1 * s3;
      dAdp[1][2][2][is] = -s1 * s2 * s3 - c1 * c3;

      dAdp[2][0][2][is] = 0;
      dAdp[2][1][2][is] = c2 * c3;
      dAdp[2][2][2][is] = -c2 * s3;

      for (Int_t i = 0; i < 3; i++)
         for (Int_t j = 0; j < 3; j++) {
            dRdp[i][j][2][is] = 0;
            for (Int_t k = 0; k < 3; k++) dRdp[i][j][2][is] += tl2g[i][k] * dAdp[k][j][2][is];
         }

      // = = = = dR=dR^(-1)/dp6dp6
      d2Adp2[0][0][2][2][is] = 0;
      d2Adp2[0][1][2][2][is] = -c1 * s2 * s3 + s1 * c3;
      d2Adp2[0][2][2][2][is] = -c1 * s2 * c3 - s1 * s3;

      d2Adp2[1][0][2][2][is] = 0;
      d2Adp2[1][1][2][2][is] = -s1 * s2 * s3 - c1 * c3;
      d2Adp2[1][2][2][2][is] = -s1 * s2 * c3 + c1 * s3;

      d2Adp2[2][0][2][2][is] = 0;
      d2Adp2[2][1][2][2][is] = -c2 * s3;
      d2Adp2[2][2][2][2][is] = -c2 * c3;

      for (Int_t i = 0; i < 3; i++)
         for (Int_t j = 0; j < 3; j++) {
            d2Rdp2[i][j][2][2][is] = 0;
            for (Int_t k = 0; k < 3; k++) d2Rdp2[i][j][2][2][is] += tl2g[i][k] * d2Adp2[k][j][2][2][is];
         }
   }

   for (Int_t i = 0; i < 3; i++)
      for (Int_t j = 0; j < 3; j++)
         if (i == j)
            dTdq[i][j] = 1;
         else
            dTdq[i][j] = 0;

   for (Int_t k = 0; k < 3; k++)
      for (Int_t i = 0; i < 6; i++)
         for (Int_t j = 0; j < 6; j++) d2Tdqdq[k][i][j] = 0;

   // print matrixes =======================================
   /*if(111==111)
   {
   for(Int_t s=0;s<12;s++) {
   printf("sector=%2d * * * * * * * R-matrixes and their derivatives * * * * * * *\n",s);
   Double_t wr[3][3];
   for(Int_t i=0;i<3;i++) for(Int_t j=0;j<3;j++) wr[i][j]=Tlo2glo[s][i][j];PrintdR("Tm1",s,wr);
   for(Int_t i=0;i<3;i++) for(Int_t j=0;j<3;j++) wr[i][j]=dRdp[i][j][0][s];PrintdR("dRd4",s,wr);
   for(Int_t i=0;i<3;i++) for(Int_t j=0;j<3;j++) wr[i][j]=d2Rdp2[i][j][0][0][s];PrintdR("d2Rd44",s,wr);
   for(Int_t i=0;i<3;i++) for(Int_t j=0;j<3;j++) wr[i][j]=d2Rdp2[i][j][0][1][s];PrintdR("d2Rd45",s,wr);
   for(Int_t i=0;i<3;i++) for(Int_t j=0;j<3;j++) wr[i][j]=d2Rdp2[i][j][1][0][s];PrintdR("d2Rd54",s,wr);
   for(Int_t i=0;i<3;i++) for(Int_t j=0;j<3;j++) wr[i][j]=d2Rdp2[i][j][0][2][s];PrintdR("d2Rd46",s,wr);
   for(Int_t i=0;i<3;i++) for(Int_t j=0;j<3;j++) wr[i][j]=d2Rdp2[i][j][2][0][s];PrintdR("d2Rd64",s,wr);
   for(Int_t i=0;i<3;i++) for(Int_t j=0;j<3;j++) wr[i][j]=dRdp[i][j][1][s];PrintdR("dRd5",s,wr);
   for(Int_t i=0;i<3;i++) for(Int_t j=0;j<3;j++) wr[i][j]=d2Rdp2[i][j][1][1][s];PrintdR("d2Rd55",s,wr);
   for(Int_t i=0;i<3;i++) for(Int_t j=0;j<3;j++) wr[i][j]=d2Rdp2[i][j][1][2][s];PrintdR("d2Rd56",s,wr);
   for(Int_t i=0;i<3;i++) for(Int_t j=0;j<3;j++) wr[i][j]=d2Rdp2[i][j][2][1][s];PrintdR("d2Rd65",s,wr);
   for(Int_t i=0;i<3;i++) for(Int_t j=0;j<3;j++) wr[i][j]=dRdp[i][j][2][s];PrintdR("dRd6",s,wr);
   for(Int_t i=0;i<3;i++) for(Int_t j=0;j<3;j++) wr[i][j]=d2Rdp2[i][j][2][2][s];PrintdR("d2Rd66",s,wr);
   }
   }*/

   // end of print    =======================================
}

//______________________________________________________________________
// set local parameters
void DerLMpdTpc::SetLocPar(Double_t *qt)
{
   // The method is called one time after the track information was read

   for (Int_t i = 0; i < 5; i++) q[i] = qt[i];
   sq4              = TMath::Sin(q[3]);
   cq4              = TMath::Cos(q[3]);
   sq5              = TMath::Sin(q[4]);
   cq5              = TMath::Cos(q[4]);
   e[0]             = sq4 * cq5;
   e[1]             = sq4 * sq5;
   e[2]             = cq4;
   dTdq[0][3]       = cq4 * cq5;
   dTdq[0][4]       = -sq4 * sq5;
   dTdq[1][3]       = cq4 * sq5;
   dTdq[1][4]       = sq4 * cq5;
   dTdq[2][3]       = -sq4;
   dTdq[2][4]       = 0;
   d2Tdqdq[0][3][3] = d2Tdqdq[0][4][4] = -sq4 * cq5;
   d2Tdqdq[0][3][4] = d2Tdqdq[0][4][3] = -cq4 * sq5;
   d2Tdqdq[1][3][3] = d2Tdqdq[1][4][4] = -sq4 * sq5;
   d2Tdqdq[1][3][4] = d2Tdqdq[1][4][3] = cq4 * cq5;
   d2Tdqdq[2][3][3]                    = -cq4;
}

//______________________________________________________________________
// set local parameters
void DerLMpdTpc::SetHitTrack(Int_t isec, Double_t *hit)
{
   // The method is called for each hit before request of partial derivatives
   // h[] - hit coordinates in GCS
   // hs[] - hit coordinates in LCS
   // T[] - coordinates of the track point closest to the hit h[] (in GCS)
   t0 = (hit[0] - q[0]) * e[0] + (hit[1] - q[1]) * e[1] + (hit[2] - q[2]) * e[2];
   for (Int_t i = 0; i < 3; i++) {
      h[i] = hit[i];
      T[i] = q[i] + e[i] * t0;
   }
   TVector3 h3(h[0], h[1], h[2]);
   TVector3 hs3 = G0_gsc[isec] + glo2loc[isec] * h3;
   hs[0]        = hs3.X();
   hs[1]        = hs3.Y();
   hs[2]        = hs3.Z();
   hs3          = G0_tlsc[isec] + Tglo2lo[isec] * h3;
   ht[0]        = hs3.X();
   ht[1]        = hs3.Y();
   ht[2]        = hs3.Z();
   h3           = TVector3(T[0], T[1], T[2]);
   hs3          = G0_tlsc[isec] + Tglo2lo[isec] * h3;
   Tt[0]        = hs3.X();
   Tt[1]        = hs3.Y();
   Tt[2]        = hs3.Z();
   // printf("hs(%f %f %f)\n",h[0],h[1],h[2]);
   // printf("ht(%f %f %f)\n",ht[0],ht[1],ht[2]);
   // printf("Tt(%f %f %f)\n",Tt[0],Tt[1],Tt[2]);
}

//______________________________________________________________________
// Get local derivatives
void DerLMpdTpc::GetLocDer(Double_t *dFdq)
{
   for (Int_t i = 0; i < 3; i++) dFdq[i] = 2 * (T[i] - h[i]);
   for (Int_t i = 3; i < 5; i++) {
      dFdq[i] = 0;
      for (Int_t n = 0; n < 3; n++) dFdq[i] += dTdq[n][i] * (T[n] - h[n]);
      dFdq[i] *= 2 * t0;
   }
}

//______________________________________________________________________
// Get 2nd order local derivatives
void DerLMpdTpc::GetLocDer2(Double_t d2Fdqdq[6][6])
{
   for (Int_t i = 0; i < 5; i++)
      for (Int_t j = 0; j < 5; j++) {
         d2Fdqdq[i][j] = 0;
         for (Int_t n = 0; n < 3; n++)
            d2Fdqdq[i][j] += 2 * (dTdq[n][i] * t0 * dTdq[n][j] * t0 + (T[n] - h[n]) * d2Tdqdq[n][i][j] * t0);
      }
}

//______________________________________________________________________
// Get global derivatives
void DerLMpdTpc::GetGloDer(Int_t isec, Double_t *dFdp)
{
   for (Int_t i = 0; i < 6; i++) {
      dFdp[i] = 0;
      for (Int_t n = 0; n < 3; n++) {
         if (i < 3) {
            dFdp[i] += Tlo2glo[isec](n, i) * (h[n] - T[n]);
            //        dFdp[i]+=(ht[n]-Tt[n]);
         } else {
            for (Int_t k = 0; k < 3; k++) {
               dFdp[i] += dRdp[n][k][i - 3][isec] * hs[k] * (h[n] - T[n]);
            }
         }
      }
      dFdp[i] *= 2;
   }
   Double_t sum[3] = {0, 0, 0};
   for (Int_t i = 0; i < 3; i++) {
      for (Int_t k = 0; k < 3; k++) sum[i] += Tlo2glo[isec](i, k) * dFdp[k];
   }
   for (Int_t i = 0; i < 3; i++) dFdp[i] = sum[i];
   // printf("  dFdp(%e %e %e %e %e %e)\n",dFdp[0],dFdp[1],dFdp[2],dFdp[3],dFdp[4],dFdp[5]);
}

/*/______________________________________________________________________
  //Get global derivatives
void DerLMpdTpc::GetGloDerT(Int_t isec,Double_t* dFdp) {
  for(Int_t i=0;i<6;i++) {
    if(i<3) dFdp[i]=(ht[i]-Tt[i]);
    else {
      dFdp[i]=0;
      for(Int_t n=0;n<3;n++) {
        for(Int_t k=0;k<3;k++) {
          dFdp[i]+=dAdp[n][k][i-3][isec]*hs[k]*(ht[n]-Tt[n]);
        }
      }
    }
    dFdp[i]*=2;
  }
//printf("dFdp_T(%e %e %e %e %e %e)\n",dFdp[0],dFdp[1],dFdp[2],dFdp[3],dFdp[4],dFdp[5]);

}*/

//______________________________________________________________________

// Get 2nd order global derivatives
void DerLMpdTpc::Getd2Fdpdp(Int_t isec, Double_t d2Fdpdp[6][6])
{
   for (Int_t i = 0; i < 6; i++) {
      for (Int_t j = i; j < 6; j++) {
         d2Fdpdp[i][j] = 0;
         if (i < 3 && j < 3) {
            for (Int_t n = 0; n < 3; n++) {
               d2Fdpdp[i][j] += Tlo2glo[isec](n, i) * Tlo2glo[isec](n, j);
            }
         } else {
            if (i < 3) {
               for (Int_t n = 0; n < 3; n++) {
                  Double_t S = 0;
                  for (Int_t k = 0; k < 3; k++) S += dRdp[n][k][j - 3][isec] * hs[k];
                  d2Fdpdp[i][j] += Tlo2glo[isec](n, i) * S;
               }
            } else {
               if (j < 3) {
                  for (Int_t n = 0; n < 3; n++) {
                     Double_t S = 0;
                     for (Int_t k = 0; k < 3; k++) S += dRdp[n][k][i - 3][isec] * hs[k];
                     d2Fdpdp[i][j] + Tlo2glo[isec](n, j) * S;
                  }
               } else {
                  for (Int_t n = 0; n < 3; n++) {
                     Double_t S1 = 0, S2 = 0, S3 = 0;
                     for (Int_t k = 0; k < 3; k++) {
                        S1 += d2Rdp2[n][k][i - 3][j - 3][isec] * hs[k];
                        S2 += dRdp[n][k][i - 3][isec] * hs[k];
                        S3 += dRdp[n][k][j - 3][isec] * hs[k];
                     }
                     d2Fdpdp[i][j] += (h[n] - T[n]) * S1 + S2 * S3;
                  }
               }
            }
         }
         d2Fdpdp[i][j] *= 2;
      }
   }
}

//______________________________________________________________________
// Get 2nd order local-global derivatives
// Before should be a call of GetGloDer(Int_t isec,Double_t* dFdp)
// to fill dHdp[6][6]
void DerLMpdTpc::GetLocGloDer2(Int_t isec, Double_t d2Fdpdq[6][6])
{
   for (Int_t i = 0; i < 5; i++) {
      Double_t t = (i > 2) ? t0 : 1;
      for (Int_t j = 0; j < 6; j++) {
         d2Fdpdq[i][j] = 0;
         for (Int_t k = 0; k < 3; k++) {
            Double_t S = 0;
            for (Int_t n = 0; n < 3; n++) S += dRdp[k][n][j - 3][isec] * hs[n];
            d2Fdpdq[i][j] += dTdq[k][i] * S * t;
            d2Fdpdq[i][j] *= -2;
         }
      }
   }
}

//______________________________________________________________________
void DerLMpdTpc::PrintdR(const char *title, Int_t s, Double_t Rm1[3][3])
{
   printf("sector=%2d =========== %s =============\n", s, title);
   for (Int_t i = 0; i < 3; i++) {
      printf("%18.10e  %18.10e  %18.10e\n", Rm1[i][0], Rm1[i][1], Rm1[i][2]);
   }
   printf("\n");
}

//______________________________________________________________________
void DerLMpdTpc::PrintM66(const char *title, Double_t m[6][6])
{
   printf("=========== %s =============\n", title);
   for (Int_t i = 0; i < 6; i++)
      printf("%18.10e  %18.10e  %18.10e  %18.10e  %18.10e  %18.10e\n", m[i][0], m[i][1], m[i][2], m[i][3], m[i][4],
             m[i][5]);
}
