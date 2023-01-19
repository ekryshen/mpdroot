void AFile(Double_t rcm, Double_t adg, Int_t Print, const char *FileName = "adir/afile.root")
{

   Double_t R0_A[3][24]; // shift  of the LSC beginning
   Double_t alpha_A[24]; // alfa angle GSC-LSC due to the alignment
   Double_t beta_A[24];  // beta angle GSC-LSC due to the alignment
   Double_t gamma_A[24]; // gamma angle GSC-LSC due to the alignment

   Double_t shift[3] = {rcm, rcm, rcm};
   Double_t grad[3]  = {adg, adg, adg};
   Double_t d2r      = TMath::Pi() / 180;
   Double_t r2d      = 180 / TMath::Pi();
   Double_t trans[3];
   for (Int_t i = 0; i < 3; i++) trans[i] = grad[i] * d2r;
   if (Print != 0) cout << "         afile: write to the file=" << FileName << endl;
   for (Int_t i = 0; i < 24; i++) {
      R0_A[0][i] = -shift[0] + 2 * shift[0] * gRandom->Rndm();
      R0_A[1][i] = -shift[1] + 2 * shift[1] * gRandom->Rndm();
      R0_A[2][i] = -shift[2] + 2 * shift[2] * gRandom->Rndm();

      alpha_A[i] = (-trans[0] + 2 * trans[0] * gRandom->Rndm());
      beta_A[i]  = (-trans[1] + 2 * trans[1] * gRandom->Rndm());
      gamma_A[i] = (-trans[2] + 2 * trans[2] * gRandom->Rndm());
      if (i < 12)
         if (Print != 0)
            printf("sector=%2d  r(%6.3f,%6.3f,%6.3f)  a(%6.3f,%6.3f,%6.3f)\n", i, R0_A[0][i], R0_A[1][i], R0_A[2][i],
                   alpha_A[i] * r2d, beta_A[i] * r2d, gamma_A[i] * r2d);
   }

   TFile *file        = new TFile(FileName, "RECREATE");
   TTree *t_alignment = new TTree("alignment", "alignment_parameters");
   t_alignment->Branch("R0_A", &R0_A, "R0_A[3][24]/D");
   t_alignment->Branch("alpha_A", &alpha_A, "alpha_A[24]/D");
   t_alignment->Branch("beta_A", &beta_A, "beta_A[24]/D");
   t_alignment->Branch("gamma_A", &gamma_A, "gamma_A[24]/D");

   t_alignment->Fill();
   t_alignment->Write();
   file->Close();
   delete file;
}
