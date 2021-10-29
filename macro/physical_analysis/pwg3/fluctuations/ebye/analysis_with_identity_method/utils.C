void tuneGraphAxisLabels( TGraph *gr )
{
    gr->GetYaxis()->SetTitleOffset( 1.45 );
    gr->GetYaxis()->SetTitleSize( 0.048 );
    gr->GetYaxis()->SetLabelSize( 0.042 );

    gr->GetXaxis()->SetTitleOffset( 0.95 );
    gr->GetXaxis()->SetTitleSize( 0.048 );
    gr->GetXaxis()->SetLabelSize( 0.042 );
}

void tuneHist1D( TH1 *hist )
{
    hist->GetYaxis()->SetTitleOffset( 1.45 );
    hist->GetYaxis()->SetTitleSize( 0.048 );
    hist->GetYaxis()->SetLabelSize( 0.042 );

    hist->GetXaxis()->SetTitleOffset( 0.95 );
    hist->GetXaxis()->SetTitleSize( 0.048 );
    hist->GetXaxis()->SetLabelSize( 0.042 );
}

void tuneHist2D( TH2 *hist )
{
    hist->GetYaxis()->SetTitleOffset( 1.45 );
    hist->GetYaxis()->SetTitleSize( 0.048 );
    hist->GetYaxis()->SetLabelSize( 0.042 );

    hist->GetXaxis()->SetTitleOffset( 0.95 );
    hist->GetXaxis()->SetTitleSize( 0.048 );
    hist->GetXaxis()->SetLabelSize( 0.042 );
}

void tuneCanvas(TCanvas *canvas)
{
    canvas->SetLeftMargin(0.15);
    canvas->SetRightMargin(0.05);
    canvas->SetTopMargin(0.05);
    canvas->SetBottomMargin(0.11);
//    canvas->SetGridy();
}

void tuneLegend( TLegend *leg )
{
    leg->SetFillColor(kWhite);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
}


void drawGraph(TGraph *gr, int marker, int color, const char* drawOpt =  "P same", double markerSize = -1, int drawCloneOrDraw = 0, double alphaFill = 1.0, double alphaLine = -1 )
{
    if ( alphaLine < 0 )
        alphaLine = alphaFill;
    gr->SetMarkerStyle( marker );
    gr->SetMarkerColorAlpha( color, alphaLine );
    gr->SetLineColorAlpha( color, alphaLine );
//    gr->SetFillColor( color );
    gr->SetFillColorAlpha(color, alphaFill);

    if (markerSize>0)
        gr->SetMarkerSize( markerSize );
    if ( drawCloneOrDraw == 0 )
        gr->DrawClone( drawOpt );
    else
        gr->Draw( drawOpt );
}

void drawHist(TH1 *hist, int marker, int color, const char* drawOpt =  "P same", double markerSize = -1, int drawCloneOrDraw = 0, double alpha = 1.0 )
{
    hist->SetMarkerStyle( marker );
    hist->SetMarkerColorAlpha( color, alpha );
    hist->SetLineColorAlpha( color, alpha );
//    hist->SetFillColor( color );
    //hist->SetFillColorAlpha(color, alpha);

    if (markerSize>0)
        hist->SetMarkerSize( markerSize );
    if ( drawCloneOrDraw == 0 )
        hist->DrawClone( drawOpt );
    else if ( drawCloneOrDraw == 1 )
        hist->Draw( drawOpt );
    else if ( drawCloneOrDraw == 2 )
        hist->DrawCopy( drawOpt );

}


void drawTex(TLatex *tex, double fontSize = 0.045, int drawCloneOrDraw = 0 ) //, int color = -1)
{
    tex->SetNDC(1);
    tex->SetTextFont(42);
    tex->SetTextSize(fontSize);
    //    tex->SetLineWidth(2);
//    if (color!=-1)
//        tex->SetTextColor( color );
    if ( drawCloneOrDraw == 0 )
        tex->DrawClone();
    else
        tex->Draw();

}

void removeErrorsForGraph( TGraphErrors *gr )
{
    for ( int i = 0; i < gr->GetN(); i++ )
        gr->SetPointError(i,0,0);
}
void calcPointsRatio(TGraphErrors *gr, TGraphErrors *grDenom, bool calcErr=false)
{
    double x, y;
    double x1, y1;
    for ( int i = 0; i < gr->GetN(); i++ )
    {
        gr->GetPoint(i,x,y);
        grDenom->GetPoint(i,x1,y1);
        double err1 = gr->GetErrorY( i );
        double err2 = grDenom->GetErrorY( i );
        double ratio = y/y1;
        gr->SetPoint(i,x,y/y1);
        if (calcErr && fabs(y)>0>0 && fabs(y1)>0)
            gr->SetPointError(i,0,sqrt(err1/y*err1/y+err2/y1*err2/y1));
    }
}

void calcPointsRatio(TGraphAsymmErrors *gr, TGraphAsymmErrors *grDenom, bool calcErr=false)
{
    double x, y;
    double x1, y1;
    for ( int i = 0; i < gr->GetN(); i++ )
    {
        gr->GetPoint(i,x,y);
        grDenom->GetPoint(i,x1,y1);
        double err1 = gr->GetErrorY( i );
        double err2 = grDenom->GetErrorY( i );
        double ratio = y/y1;
        gr->SetPoint(i,x,y/y1);
        if (calcErr && fabs(y)>0>0 && fabs(y1)>0)
        {
            double newErr = sqrt(err1/y*err1/y+err2/y1*err2/y1);
            gr->SetPointError(i,0,0,newErr,newErr);
        }
    }
}


void multiplyPointsByFactor(TGraphErrors *gr, double factor, bool flagMultErrorY=false )
{
    double x, y;
    for ( int i = 0; i < gr->GetN(); i++ )
    {
        gr->GetPoint(i,x,y);
        gr->SetPoint(i,x,y*factor);
        if ( flagMultErrorY )
            gr->SetPointError(i,gr->GetErrorX(i),gr->GetErrorY(i)*factor);
    }
}

void squarePoints(TGraphErrors *gr, bool flagMultErrorY=false )
{
    double x, y;
    for ( int i = 0; i < gr->GetN(); i++ )
    {
        gr->GetPoint(i,x,y);
        gr->SetPoint(i,x,y*y);
        if ( flagMultErrorY )
            gr->SetPointError( i, gr->GetErrorX(i),gr->GetErrorY(i)*gr->GetErrorY(i) );
    }
}

void multiplyPointsByFactorsFromAnotherGraph(TGraphErrors *gr, TGraphErrors *grFactors ) //, bool flagMultErrorY=false )
{
    double x, y;
    double x1, factor;
    for ( int i = 0; i < gr->GetN(); i++ )
    {
        gr->GetPoint(i,x,y);
        grFactors->GetPoint(i,x1,factor);

        gr->SetPoint(i,x,y*factor);
//        if ( flagMultErrorY )
            gr->SetPointError(i,gr->GetErrorX(i),gr->GetErrorY(i)*1.0/*factor*/);
    }
}

void addToPoints(TGraphErrors *gr, double add)
{
    double x, y;
    for ( int i = 0; i < gr->GetN(); i++ )
    {
        gr->GetPoint(i,x,y);
        gr->SetPoint(i,x,y+add);
    }
}


void reflectAndCombinePoints( TGraphErrors *gr )
{
    double x1, y1, x2, y2;
    int nP = gr->GetN();
    for ( int i = 1; i < nP/2; i++ )
    {
        int reflPoint = nP-1 - (i-1);
        gr->GetPoint( i, x1, y1 );
        gr->GetPoint( reflPoint, x2, y2 );
        gr->SetPoint(i, x1, (y1+y2)/2 ) ;

        // errors:
        double err1 = gr->GetErrorY( i );
        double err2 = gr->GetErrorY( reflPoint );

        gr->SetPointError( i, 0, (err1+err2)/2. );
    }
}



void getQuantiles(TH1 *h, const int nq, double *yq, bool printOrNot = true)
{
    //    const Int_t nq = 10;
    //    const Int_t nshots = 10;
    Double_t xq[nq];  // position where to compute the quantiles in [0,1]
    //    Double_t yq[nq];  // array to contain the quantiles
    for (Int_t i=0;i<nq;i++)
        xq[i] = Float_t(i+1)/nq;

    h->GetQuantiles( nq, yq, xq );

    if (printOrNot)
    {
        cout << "nq=" << nq << endl;
        for (Int_t i=0;i<nq;i++)
            cout << yq[i] << ", ";
        cout << endl;
    }
}

void rearrangeBoundaries( const int nq, double *multBounds, double *boundsMin, double *boundsMax)
{
    for (Int_t i=0;i<nq;i++)
    {
        if ( i == -1 )
            boundsMin[i] = 0;
        else
            boundsMin[i] = multBounds[i-1];
        boundsMax[i] = multBounds[i+1];
        cout << "bounds for i=" << i << ": " << boundsMin[i] << " " << boundsMax[i] << endl;
    }
}


void drawCanvasWithClasses(TH1D *hist1D, /*TString label,*/ const int nCentrBins, double *centrMultBounds, double *multBinCenters)//, TCanvas *canv_mult_withClasses )
{
//    TCanvas *
//    canv_mult_withClasses = new TCanvas(Form("canv_mult_withClasses_%s",label.Data())
//                                                 ,Form("canv_mult_withClasses_%s",label.Data()),300,100,800,600 );
//    tuneCanvas(canv_mult_withClasses);

    hist1D->SetLineColor(kRed);
    hist1D->DrawCopy();

    gPad->SetLogy();

    //draw centrality classes on mult hist
    //    int nBinsForClassHists = fHistParticlesInCutConditionInEvent->GetNbinsX();
    TH1D **fHistCentrClass = new TH1D*[nCentrBins];
    for ( int iCentrClass = 0; iCentrClass < nCentrBins; iCentrClass++ )
    {
        fHistCentrClass[iCentrClass] = (TH1D*)hist1D->Clone( Form("fHistCentrClass%d", iCentrClass) );   //new TH1D( Form("fHistCentrClass%d", iCentrClass), Form("iCentrClass%d", iCentrClass)
        //                                                 , hist1D->GetNbinsX(), hist1D->GetBinLowEdge(1), hist1D->GetBin );
        fHistCentrClass[iCentrClass]->Reset();
    }

    int centrClassId = 0;
    for ( int iBin = 0; iBin < hist1D->GetNbinsX(); iBin++ )
    {
        //        cout << fHistParticlesInCutConditionInEvent->GetBinCenter(iBin+1) << " " << centrMultBounds[centrClassId] << endl;
        if ( hist1D->GetBinCenter(iBin+1) > centrMultBounds[centrClassId] )
            centrClassId++;
        if ( centrClassId >= nCentrBins )
            break;
        double binContent = hist1D->GetBinContent( iBin+1 );
        fHistCentrClass[centrClassId]->SetBinContent( iBin+1, binContent );
    }
    for ( int iCentrClass = 0; iCentrClass < nCentrBins; iCentrClass++ )
    {
        fHistCentrClass[iCentrClass]->SetFillColor( kOrange - 5 + iCentrClass );
        fHistCentrClass[iCentrClass]->SetLineColor(kBlue);
        fHistCentrClass[iCentrClass]->DrawCopy( "same" );
        multBinCenters[iCentrClass] = fHistCentrClass[iCentrClass]->GetMean();

        //                double meanNchInCentrBin = fHistCentrClass[iCentrClass]->GetMean();
        //                fHistMultClassMeanNch->SetBinContent( iCentrClass+1, meanNchInCentrBin );
        //                cout << "meanNch in centrality bin " << iCentrClass << ": " << meanNchInCentrBin << endl;
    }

    //canv_mult_withClasses->SaveAs( Form("output/canv_%s_classes_%d.eps", label.Data(), nCentrBins));

}


void multiplyXbyFactor(TGraphErrors *gr, double factor )
{
    double x, y;
    for ( int i = 0; i < gr->GetN(); i++ )
    {
        gr->GetPoint(i,x,y);
        gr->SetPoint(i,x*factor,y);
    }
}

void shiftPointX(TGraphErrors *gr, double shift)
{
    double x, y;
    for ( int i = 0; i < gr->GetN(); i++ )
    {
        gr->GetPoint(i,x,y);
        gr->SetPoint(i,x+shift,y);
    }
}
