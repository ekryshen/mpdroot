#include <Rtypes.h>

#if !defined(__CINT__) && !defined(__CLING__)
// ROOT includes
#include <TFitter.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TGraph.h>
#include <TGraph2D.h>
#include <TProfile.h>
#include <TString.h>
#include <TLatex.h>
#include <TFile.h>
#include <TLeaf.h>
#include <TTree.h>
#include <TNtuple.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TStopwatch.h>
#include <TSystem.h>
#include <TChain.h>

using namespace std;

#include <stdlib.h>
#include <stdio.h>
#include <iostream>

//void mntfnc(int&, double*, double&, double*, int);
TCanvas* plot(const char*, char*);

int main()
{
    if (!plot("../output.root", "plots.root")) return 1;
    return 0;
}

#else // __CINT__

R__ADD_INCLUDE_PATH($VMCWORKDIR)
#include "macro/mpd/mpdloadlibs.C"

#endif

//#define _DEBUG

const unsigned nSectors = 24;


struct OutStruct
{
    Double_t t[nSectors];
    Double_t ann[nSectors];
    Double_t nsg[nSectors];
    Double_t nht[nSectors];
    Double_t dx[nSectors];
    Double_t dy[nSectors];
    Double_t dz[nSectors];
    Double_t dp[nSectors];
    Double_t de[nSectors];
    Double_t dk[nSectors];
    Double_t Dx[nSectors];
    Double_t Dy[nSectors];
    Double_t Dz[nSectors];
    Double_t Dp[nSectors];
    Double_t De[nSectors];
    Double_t Dk[nSectors];
};


TCanvas* plot(const char* fname_inp = "../output.root", char* fname_out = nullptr)
{
    auto *file_inp = TFile::Open(fname_inp);

    if (!file_inp)
    {
        cout << "Error opening input data: " << fname_inp << endl;
        return nullptr;
    }

    auto *tree_inp = (TTree *) file_inp->Get("tree");

    if (!tree_inp)
    {
        cout << "Null pointer to input tree structure retreived" << endl;
        return nullptr;
    }

    OutStruct st;
    TBranch *br;

    tree_inp->SetBranchAddress("st", st.t, &br);
    auto n_ent = tree_inp->GetEntriesFast();

    cout << "Events: " << n_ent << endl;

    Double_t *step = new Double_t[n_ent];
    std::vector<Double_t>  dx_tot, sx_tot, errx_tot, dy_tot, sy_tot, erry_tot, dz_tot, sz_tot, errz_tot;
    std::vector<Double_t>  ann[nSectors], dx[nSectors], sx[nSectors], errx[nSectors],
                                          dy[nSectors], sy[nSectors], erry[nSectors],
                                          dz[nSectors], sz[nSectors], errz[nSectors];
    Double_t ind[nSectors], x_0[nSectors], dx_min[nSectors], sx_min[nSectors], errx_min[nSectors],
                            y_0[nSectors], dy_min[nSectors], sy_min[nSectors], erry_min[nSectors],
                            z_0[nSectors], dz_min[nSectors], sz_min[nSectors], errz_min[nSectors];

    for (unsigned k = 0; k < nSectors; ++k)
        dx_min[k] = errx_min[k] = dy_min[k] = erry_min[k] = dz_min[k] = errz_min[k] = 1.E10;

    auto prec = 0.3;
    auto sec_select = 1;//1; // each sec_select sector will be shown

    for (auto i = 0; i < n_ent; ++i)
    {
        step[i] = i;
        tree_inp->GetEntry(i);

        for (unsigned k = 0; k < nSectors; ++k)
        {
            ann[k].push_back(st.ann[k]);

            if (!i)
                x_0[k] = -st.dx[k];

            auto diff = st.dx[k];
            auto err  = fabs(diff)/fabs(x_0[k]);
            auto val  = st.Dx[k];

            dx[k].push_back(diff);
            sx[k].push_back(val);
            errx[k].push_back(err);

            if (err < prec)
            {
                dx_tot.push_back(diff);
                sx_tot.push_back(val);
                errx_tot.push_back(err);
            }

            if (fabs(diff) < fabs(dx_min[k]))
            {
                dx_min[k] = diff;
                sx_min[k] = val;
            }

            if (!i)
                y_0[k] = -st.dy[k];

            diff = st.dy[k];
            err  = fabs(diff)/fabs(y_0[k]);
            val  = st.Dy[k];

            dy[k].push_back(diff);
            sy[k].push_back(val);
            erry[k].push_back(err);

            if (err < prec)
            {
                dy_tot.push_back(diff);
                sy_tot.push_back(val);
                erry_tot.push_back(err);
            }

            if (fabs(diff) < fabs(dy_min[k]))
            {
                dy_min[k] = diff;
                sy_min[k] = val;
            }

            if (!i)
                z_0[k] = -st.dz[k];

            diff = st.dz[k];
            err  = fabs(diff)/fabs(z_0[k]);
            val  = st.Dz[k];

            dz[k].push_back(diff);
            sz[k].push_back(val);
            errz[k].push_back(err);

            if (err < prec)
            {
                dz_tot.push_back(diff);
                sz_tot.push_back(val);
                errz_tot.push_back(err);
            }

            if (fabs(diff) < fabs(dz_min[k]))
            {
                dz_min[k] = diff;
                sz_min[k] = val;
            }
        }
    }

    for (unsigned k = 0; k < nSectors; ++k)
    {
        ind[k] = k;

        auto diff = fabs(dx_min[k]);
        auto err  = diff/fabs(x_0[k]);
        auto val  = sx_min[k];

        //if (err > prec)
        //    cout << "sector " << k << " X BAD CONV: min diff " << diff << " Dx " <<  val << endl;
        errx_min[k] = err;

        diff = fabs(dy_min[k]);
        err  = diff/fabs(y_0[k]);
        val  = sy_min[k];

        //if (err > prec)
        //    cout << "sector " << k << " Y BAD CONV: min diff " << diff << " Dy " <<  val << endl;
        erry_min[k] = err;

        diff = fabs(dz_min[k]);
        err  = diff/fabs(z_0[k]);
        val  = sz_min[k];

        //if (err > prec)
        //    cout << "sector " << k << " Z BAD CONV: min diff " << diff << " Dz " <<  val << endl;
        errz_min[k] = err;
    }

    char buf[255];

    if (!fname_out)
    {
        fname_out = new char[255];
        auto num = strlen(fname_inp)-5;

        sprintf(buf, "_plots.root");
        memmove(fname_out, fname_inp, num);
        memmove(fname_out+num, buf, strlen(buf)+1);
    }

    TFile *file_out = TFile::Open(fname_out, "RECREATE");

    TGraph *grxv[nSectors], *gryv[nSectors], *grzv[nSectors], *grx[nSectors], *gry[nSectors], *grz[nSectors];//  *grdsx[nSectors],  *grdsy[nSectors];
    TCanvas *cx[nSectors], *cy[nSectors], *cz[nSectors];//, *cds[nSectors];

    for (unsigned k = 0; k < nSectors; ++k)
    {
        if (k%sec_select) continue;

        // X diff

        sprintf(buf, "x_diff_vs_cov_%02d", k);
        cx[k] = new TCanvas(buf, buf, 0, 0, 1600, 800);
        cx[k]->Divide(2,1);

        grx[k] = new TGraph(sx[k].size(), step, sx[k].begin().base());

        cx[k]->cd(2)->SetGridx();
        grx[k]->SetMarkerStyle(20);
        grx[k]->SetMarkerSize(0.4);
        //grx[k]->Draw("A RX P");
        grx[k]->Draw("AP");
        cx[k]->Update();

        grx[k]->SetTitle(buf);
        grx[k]->GetXaxis()->SetTitle("step");
        grx[k]->GetYaxis()->SetTitle("Dx");

        sprintf(buf, "x_diff_vs_step_%02d", k);
        grxv[k] = new TGraph(dx[k].size(), step, dx[k].begin().base());

        cx[k]->cd(1)->SetGridx();
        grxv[k]->SetMarkerStyle(20);
        grxv[k]->SetMarkerSize(0.4);
        grxv[k]->Draw("AP");
        cx[k]->Update();

        grxv[k]->SetTitle(buf);
        grxv[k]->GetYaxis()->SetTitle("x diff");
        grxv[k]->GetXaxis()->SetTitle("step");

        //sprintf(buf, "x_diff_vs_cov_%02d.png", k);
        //cx[k]->Print(buf);

        // Y diff

        sprintf(buf, "y_diff_vs_cov_%02d", k);
        cy[k] = new TCanvas(buf, buf, 0, 0, 1600, 800);
        cy[k]->Divide(2,1);

        gry[k] = new TGraph(sy[k].size(), step, sy[k].begin().base());

        cy[k]->cd(2)->SetGridx();
        gry[k]->SetMarkerStyle(20);
        gry[k]->SetMarkerSize(0.4);
        gry[k]->Draw("AP");
        cy[k]->Update();

        gry[k]->SetTitle(buf);
        gry[k]->GetXaxis()->SetTitle("step");
        gry[k]->GetYaxis()->SetTitle("Dy");

        sprintf(buf, "y_diff_vs_step_%02d", k);
        gryv[k] = new TGraph(dy[k].size(), step, dy[k].begin().base());

        cy[k]->cd(1)->SetGridx();
        gryv[k]->SetMarkerStyle(20);
        gryv[k]->SetMarkerSize(0.4);
        gryv[k]->Draw("AP");
        cy[k]->Update();

        gryv[k]->SetTitle(buf);
        gryv[k]->GetYaxis()->SetTitle("y diff");
        gryv[k]->GetXaxis()->SetTitle("step");

        //sprintf(buf, "y_diff_vs_cov_%02d.png", k);
        //cy[k]->Print(buf);

        // Z diff

        sprintf(buf, "z_diff_vs_cov_%02d", k);
        cz[k] = new TCanvas(buf, buf, 0, 0, 1600, 800);
        cz[k]->Divide(2,1);

        grz[k] = new TGraph(sz[k].size(), step, sz[k].begin().base());

        cz[k]->cd(2)->SetGridx();
        grz[k]->SetMarkerStyle(20);
        grz[k]->SetMarkerSize(0.4);
        grz[k]->Draw("AP");
        cz[k]->Update();

        grz[k]->SetTitle(buf);
        grz[k]->GetXaxis()->SetTitle("step");
        grz[k]->GetYaxis()->SetTitle("Dz");

        sprintf(buf, "z_diff_vs_step_%02d", k);
        grzv[k] = new TGraph(dz[k].size(), step, dz[k].begin().base());

        cz[k]->cd(1)->SetGridx();
        grzv[k]->SetMarkerStyle(20);
        grzv[k]->SetMarkerSize(0.4);
        grzv[k]->Draw("AP");
        cz[k]->Update();

        grzv[k]->SetTitle(buf);
        grzv[k]->GetYaxis()->SetTitle("z diff");
        grzv[k]->GetXaxis()->SetTitle("step");

        //sprintf(buf, "z_diff_vs_cov_%02d.png", k);
        //cz[k]->Print(buf);

        /*
        sprintf(buf, "x_cov_vs_step_%02d", k);
        cds[k] = new TCanvas(buf, buf, 0, 0, 1600, 800);
        cds[k]->Divide(2,1);

        grdsx[k] = new TGraph(sx[k].size(), step, sx[k].begin().base());

        cds[k]->cd(1)->SetGridx();
        cds[k]->cd(1)->SetGridy();
        cds[k]->cd(1)->SetLogy();
        grdsx[k]->SetMarkerStyle(20);
        grdsx[k]->SetMarkerSize(0.4);
        grdsx[k]->Draw("AP");
        cds[k]->Update();

        grdsx[k]->SetTitle(buf);
        grdsx[k]->GetYaxis()->SetTitle("Dx");
        grdsx[k]->GetXaxis()->SetTitle("step");

        sprintf(buf, "y_cov_vs_step_%02d", k);
        grdsy[k] = new TGraph(sy[k].size(), step, sy[k].begin().base());

        cds[k]->cd(2)->SetGridx();
        cds[k]->cd(2)->SetGridy();
        cds[k]->cd(2)->SetLogy();
        grdsy[k]->SetMarkerStyle(20);
        grdsy[k]->SetMarkerSize(0.4);
        grdsy[k]->Draw("AP");
        cds[k]->Update();

        grdsy[k]->SetTitle(buf);
        grdsy[k]->GetYaxis()->SetTitle("Dy");
        grdsy[k]->GetXaxis()->SetTitle("step");
        */
        //sprintf(buf, "cov_vs_step_%02d.png", k);
        //cds[k]->Print(buf);
    }

    /*
    sprintf(buf, "diff_vs_cov");
    auto* c0 = new TCanvas(buf, buf, 0, 0, 1600, 800);
    c0->Divide(2, 1);

    auto* grxt  = new TGraph(dx_tot.size(), sx_tot.begin().base(), dx_tot.begin().base());
    auto* gryt  = new TGraph(dy_tot.size(), sy_tot.begin().base(), dy_tot.begin().base());

    c0->cd(1)->SetGridx();
    //c0->cd(1)->SetLogx();
    grxt->SetMarkerStyle(20);
    grxt->SetMarkerSize(0.4);
    grxt->Draw("A RX P");
    c0->Update();

    grxt->SetTitle("X");
    grxt->GetYaxis()->SetTitle("x diff");
    grxt->GetXaxis()->SetTitle("Dx");

    c0->cd(2)->SetGridx();
    //c0->cd(2)->SetLogx();
    gryt->SetMarkerStyle(20);
    gryt->SetMarkerSize(0.4);
    gryt->Draw("A RX P");
    c0->Update();

    gryt->SetTitle("Y");
    gryt->GetYaxis()->SetTitle("y diff");
    gryt->GetXaxis()->SetTitle("Dy");

    sprintf(buf, "diff_vs_cov.png");
    //c0->Print(buf);

    sprintf(buf, "err_vs_cov");
    auto* c3 = new TCanvas(buf, buf, 0, 0, 1600, 800);
    c3->Divide(2, 1);

    auto* grxet  = new TGraph(errx_tot.size(), sx_tot.begin().base(), errx_tot.begin().base());
    auto* gryet  = new TGraph(erry_tot.size(), sy_tot.begin().base(), erry_tot.begin().base());

    c3->cd(1)->SetGridx();
    grxet->SetMarkerStyle(20);
    grxet->SetMarkerSize(0.4);
    grxet->Draw("A RX P");
    c3->Update();

    grxet->SetTitle("X");
    grxet->GetYaxis()->SetTitle("x err");
    grxet->GetXaxis()->SetTitle("Dx");

    c3->cd(2)->SetGridx();
    gryet->SetMarkerStyle(20);
    gryet->SetMarkerSize(0.4);
    gryet->Draw("A RX P");
    c3->Update();

    gryet->SetTitle("Y");
    gryet->GetYaxis()->SetTitle("y err");
    gryet->GetXaxis()->SetTitle("Dy");

    sprintf(buf, "err_vs_cov.png");
    //c3->Print(buf);

    sprintf(buf, "min_err_and_diff");
    auto* c1 = new TCanvas(buf, buf, 0, 0, 1600, 800);
    c1->Divide(2, 2);

    auto* grxe = new TGraph(nSectors, ind, errx_min);
    auto* grye = new TGraph(nSectors, ind, erry_min);

    c1->cd(1);
    grxe->SetMarkerStyle(34);
    grxe->SetMarkerSize(1.5);
    grxe->Draw("AP");
    c1->Update();

    grxe->SetTitle("X");
    grxe->GetXaxis()->SetTitle("sector");
    grxe->GetYaxis()->SetTitle("x err");

    c1->cd(2);
    grye->SetMarkerStyle(34);
    grye->SetMarkerSize(1.5);
    grye->Draw("AP");
    c1->Update();

    grye->SetTitle("Y");
    grye->GetXaxis()->SetTitle("sector");
    grye->GetYaxis()->SetTitle("y err");

    auto* grxm = new TGraph(nSectors, sx_min, dx_min);
    auto* grym = new TGraph(nSectors, sy_min, dy_min);

    c1->cd(3);
    grxm->SetMarkerStyle(20);
    grxm->SetMarkerSize(0.6);
    grxm->Draw("A RX P");
    c1->Update();

    grxm->SetTitle("X");
    grxm->GetYaxis()->SetTitle("x diff min");
    grxm->GetXaxis()->SetTitle("Dx");

    c1->cd(4);
    grym->SetMarkerStyle(20);
    grym->SetMarkerSize(0.6);
    grym->Draw("A RX P");
    c1->Update();

    grym->SetTitle("Y");
    grym->GetYaxis()->SetTitle("y diff min");
    grym->GetXaxis()->SetTitle("Dy");

    sprintf(buf, "min_err_and_diff.png");
    //c1->Print(buf);

    sprintf(buf, "min_covy_vs_covx");
    auto* cc = new TCanvas(buf, buf, 0, 0, 1600, 800);
    auto* grcc = new TGraph(nSectors, sx_min, sy_min);

    cc->cd();
    grcc->SetMarkerStyle(20);
    grcc->SetMarkerSize(0.6);
    grcc->Draw("AP");
    cc->Update();

    grcc->SetTitle("min cov diag xy");
    grcc->GetXaxis()->SetTitle("Dx");
    grcc->GetYaxis()->SetTitle("Dy");

    sprintf(buf, "min_covy_vs_covx.png");
    //cc->Print(buf);

    sprintf(buf, "min_diffy_vs_diffx");
    auto* cd = new TCanvas(buf, buf, 0, 0, 1600, 800);
    auto* grdd = new TGraph(nSectors, dx_min, dy_min);

    cd->cd();
    grdd->SetMarkerStyle(20);
    grdd->SetMarkerSize(0.6);
    grdd->Draw("AP");
    cd->Update();

    grdd->SetTitle("min diff xy");
    grdd->GetXaxis()->SetTitle("diff x");
    grdd->GetYaxis()->SetTitle("diff y");

    sprintf(buf, "min_diffy_vs_diffx.png");
    //cd->Print(buf);

    sprintf(buf, "min_erry_vs_errx");
    auto* ce = new TCanvas(buf, buf, 0, 0, 1600, 800);
    auto* gree = new TGraph(nSectors, errx_min, erry_min);

    ce->cd();
    gree->SetMarkerStyle(20);
    gree->SetMarkerSize(0.6);
    gree->Draw("AP");
    ce->Update();

    gree->SetTitle("min err xy");
    gree->GetXaxis()->SetTitle("rel err x");
    gree->GetYaxis()->SetTitle("rel err y");

    sprintf(buf, "min_erry_vs_errx.png");
    //ce->Print(buf);

    c0->Write();
    c1->Write();
    c3->Write();
    cc->Write();
    cd->Write();
    ce->Write();
*/

    for (unsigned k = 0; k < nSectors; ++k)
    {
        if (k%sec_select) continue;

        cx[k]->Write(); cy[k]->Write(); cz[k]->Write();
        //cds[k]->Write();
    }

    file_out->Write();
    file_out->Close();

#ifdef DEBUG

    TH1D *hh1 = new TH1D("hh1", "hit-to-track dz", 100, (Axis_t) -1., (Axis_t) 1.);
    TH1D *hh2 = new TH1D("hh2", "hit-to-track dr (xy plane)", 100, (Axis_t) -1., (Axis_t) 1.);
    TH1D *hh3 = new TH1D("hh3", "hits R", 100, (Axis_t) 35., (Axis_t) 105.);
    TH2D *hh4 = new TH2D("hh4", "track hits dr vs R (xy plane)", 100, (Axis_t) -1., (Axis_t) 1., 100, (Axis_t) 35., (Axis_t) 105.);
    TH2D *hh5 = new TH2D("hh5", "xy planes", 100, (Axis_t) -125., (Axis_t) 125., 100, (Axis_t) -125., (Axis_t) 125.);

    c0 = new TCanvas("c0","debug", 0, 0, 1200, 1200);
    c0->Divide(2,2);

    c0->cd(1);
    hh4->GetXaxis()->SetTitle("cm");
    hh4->GetYaxis()->SetTitle("cm");
    hh4->Draw("cont");

    TGraph2D *gr = new TGraph2D(2);
    gr->SetPoint(0, -130., -130., -170.);
    gr->SetPoint(1,  130.,  130.,  170.);
    gr->SetTitle("tracks; X, cm; Y, cm; Z, cm");
    gr->Draw("P");

    for (Int_t k = 0; k < trks.size(); ++k)
    {
        auto n = trks[k].n;

        if (n)
        {
            TGraph2D *graph = new TGraph2D(n);

            for (Int_t l = 0; l < n; ++l)
                graph->SetPoint(l, trks[k].x[l], trks[k].y[l], trks[k].z[l]);

            graph->Draw("same LINE");
        }
    }

#endif

    return nullptr;
}
