//Preliminary graphs of quantities
//
//
// c++ -o pgraphs pgraphs.cpp `root-config --cflags --glibs`
//./pt HH.root ttbar.root

#include <iostream>
#include "TH1.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TTree.h"
#include <TStyle.h>
#include <TMath.h>
#include <TF1.h>
#include <TLegend.h>
#include <THStack.h>
#include <TApplication.h>
#include <fstream>

using namespace std;

//kolmogrov test
//minore è il valore di alpha più le cdf si discostano.
//Alpha rappresenta la probabilità di riscontrare una distanza (dist) supponsendo che i due campioni di dati seguano la stessa distribuzione
//Size è supposta uguale per entrambe le distribuzioni (pari al numero di bin)

void kolmogrov(TH1F * h1, TH1F * h2,  int size){
	long double alpha = 0.0001, dist = 0 ;
	for(int j = 1; j <= size; j++){
		double v1 = h1->GetBinContent(j);
		double v2 = h2->GetBinContent(j);
		if( abs(v1-v2) >= dist ) dist = abs(v1-v2);
	}
	//long double c = sqrt(-0.5*log(alpha/2.));
	//cout << "Max distance: "<<dist << " Test limit:" << c*pow(2./size, 0.5) << endl;
	alpha = 2. * exp(-1.*pow(dist ,2.)*size);	
	cout<< "The null hypostesis is rejected with a confidence level of:" << alpha *100. << "%" << endl;
	return;
}

void setstack(TH1F * hs, TH1F * hb, THStack * stack ){	
	hb->SetFillStyle(3003);
	hs->SetFillStyle(3004);
    hs->SetFillColor(2);
	hb->SetFillColor(4);
	hs->SetLineColor(2);
	hb->SetLineColor(4);	
	stack->Add(hs, "S");
	stack->Add(hb, "S");
    stack->Draw("nostack");
	gPad->SetGrid(1,1);		
	TLegend *legend1 = new TLegend(0.8,0.2,0.98,0.38);
    legend1->AddEntry(hs,"Signal", "f");
    legend1->AddEntry(hb,"Background", "f");
    legend1->Draw("SAME");
	return;
}

void DrawHisto (TTree * signal, TTree * background, string * var, TCanvas ** c, TH1F *** h, int N ){	

	for(int i = 0; i < N; i++){
		//INSERIRE RANGE DEGLI ISTOGRAMMI
		double min, max;
		int nbin;
		cout << "Inserire range della grandezza: " << var[i] <<endl;
		cout << "Min: "; cin >> min;
		cout << "Max: "; cin >> max;
		cout << "Nbin: "; cin >> nbin; 

		string title = "Graph_" + var[i];
		h[i] = new TH1F * [13];
		c[i] = new TCanvas(Form("c%d",i), title.c_str() );	
		c[i]->Divide(3,3);

		//SCRITTURA ISTOGRAMMA DISTRIBUZIONE
		c[i]->cd(1);
		string title1 = "Background_" + var[i];
		string title2 = "Signal_" + var[i];
		string var1 = var[i] + ">>" + title1;
		string var2 = var[i] + ">>" + title2;
		THStack * dist = new THStack("dist", "Distribution histo");
   		h[i][1] = new TH1F( title2.c_str(), "Signal", nbin,min,max);
 		h[i][2] = new TH1F( title1.c_str(), "Background", nbin,min,max); 
   		signal->Draw( var2.c_str(), "" );
    	background->Draw( var1.c_str(), "" );
		setstack( h[i][1], h[i][2], dist );
    	dist->SetTitle("Distribution; Energy [GeV] ; Events");

   		//ISTOGRAMMA pdf
		c[i]->cd(2);
		THStack * pdf = new THStack("pdf", "Pdf histo");
   		h[i][3] = new TH1F(Form("ps%d", i),"Signal Probability distribution", nbin,min,max);
		h[i][4] = new TH1F(Form("pb%d", i),"Background Probability distribution", nbin,min,max);
     	for(int j = 1; j <= nbin; j++){
    		h[i][3]->SetBinContent(j, h[i][1]->GetBinContent(j)/h[i][1]->GetEntries());
			h[i][4]->SetBinContent(j, h[i][2]->GetBinContent(j)/h[i][2]->GetEntries());
		}
		
		setstack( h[i][3], h[i][4], pdf );
		pdf->SetTitle("Probability distribution; Energy [GeV] ; Prob");
//		cout<<"Total probability signal: "<< h[i][3]->Integral() <<"\n"; 
//		cout<<"Total probability background: "<< h[i][4]->Integral() <<"\n"; 	
	
		//ISTOGRAMMA cdf CALCOLATA A PARTIRE DA ESTREMO SX
		c[i]->cd(3);
		THStack * cdf = new THStack("cdf", "Cdf histo");
		h[i][5]= new TH1F(Form("cps%d", i),"Signal Cumulative Probability distribution", nbin,min,max);
		h[i][6] = new TH1F(Form("cpb%d", i),"Background Cumulative Probability distribution", nbin,min,max);
    	for(int j = 1; j <= nbin; j++){
    		h[i][5]->SetBinContent(j, h[i][3]->Integral(1,j,""));
			h[i][6]->SetBinContent(j, h[i][4]->Integral(1,j,""));
		}
		setstack( h[i][5], h[i][6], cdf );
		cdf->SetTitle("Cumulative Probability distribution; Energy [GeV] ; Prob");

		//FATTORI DI MERITO CALCOLATI DA SX
    	h[i][7] = new TH1F(Form("S-B %d", i)," T factor B-S sx", nbin, min, max);
		h[i][8] = new TH1F(Form("S/B %d", i),"T factor B/S sx", nbin, min, max);
		h[i][9] = new TH1F(Form("S/sqrt(B) %d", i),"T factor B/sqrt(S) sx ", nbin, min, max);
    	for(int j = 1; j <= nbin; j++){
			double s_int = h[i][5]->GetBinContent(j);
			double b_int = h[i][6]->GetBinContent(j);
    		h[i][7]->SetBinContent(j, s_int - b_int);
			h[i][8]->SetBinContent(j, s_int / b_int);
			h[i][9]->SetBinContent(j, s_int / sqrt(b_int) );
		}
		kolmogrov( h[i][5], h[i][6], nbin);
		for (int j = 0; j < 3; j++){
			c[i]->cd(4+j);
			gPad->SetGrid(1,1);
			gStyle->SetOptStat(0000);	
			h[i][7+j]->SetLineWidth(3);
			h[i][7+j]->SetLineColor(2);
			h[i][7+j]->Draw("C");
			if (j == 0 ) h[i][7]->SetTitle("T factor S-B sx; cut; T ");
    		if (j == 1 ) h[i][8]->SetTitle("T factor S/B sx; cut; T ");
			if (j == 2 ) h[i][9]->SetTitle("T factor S/sqrt(B) sx; cut; T ");	
		}

		//FATTORI DI MERITO CALCOLATI DA dX
    	h[i][10] = new TH1F(Form("S-B dx %d", i)," T factor B-S", nbin, min, max);
		h[i][11] = new TH1F(Form("S/B  dx %d", i),"T factor B/S", nbin, min, max);
		h[i][12] = new TH1F(Form("S/sqrt(B) dx %d", i),"T factor B/sqrt(S)", nbin, min, max);
    	for(int j = 1; j <= nbin; j++){
			double s_dx = h[i][3]->Integral(j,nbin,"");
			double b_dx = h[i][4]->Integral(j,nbin,"");
    		h[i][10]->SetBinContent(j, s_dx - b_dx);
			h[i][11]->SetBinContent(j, s_dx / b_dx);
			h[i][12]->SetBinContent(j, s_dx / sqrt(b_dx) );
		}
		for (int j = 0; j < 3; j++){
			c[i]->cd(7+j);
			gPad->SetGrid(1,1);
			gStyle->SetOptStat(0000);	
			h[i][10+j]->SetLineWidth(3);
			h[i][10+j]->SetLineColor(2);
			h[i][10+j]->Draw("C");
			if (j == 0 ) h[i][10]->SetTitle("T factor S-B dx; cut; T ");
    		if (j == 1 ) h[i][11]->SetTitle("T factor S/B dx; cut; T ");
			if (j == 2 ) h[i][12]->SetTitle("T factor S/sqrt(B) dx; cut; T ");	
		}

	}
	return;
}

//Le grandezze di cui fare i grafici vanno passate sulla riga di comando.
//Per vedere i nomi con cui sono salvate nei TTree, consultare il file di lettura reading.cpp
//Il programma chiede di inserire il range di ogni variabile. Le grandezze hanno  circa i seguenti range:
//	var_phi: -pi, pi
//	var_eta: -4, 4
//	var_Et: 0, 350
//	var_pt:	0, 350
//	Ht, Htnu: 120, 800
//	Ptm: 0, 150

int main(int argc, char** argv){
    if (argc < 1){
        cout << "Usage: " << argv[0] << " q1 q2 q3 ... " << endl;
        return 1;
    } s

	TApplication * Grafica = new TApplication("App", 0, 0);
	//LETTURA DEI TTree
    TFile * sinput = TFile::Open( "HH.root" );
    TFile * binput = TFile::Open( "ttbar.root" );
    TTree * signal  = (TTree*)sinput->Get("tree");
    TTree * background  = (TTree*)binput->Get("tree");
	
	int N = argc-1;
	TCanvas ** c = new TCanvas * [N];
	TH1F *** h = new  TH1F ** [N];
	string * var = new string [N]; 
	for( int i = 0; i < N; i++) var[i] = argv[i+1];
	DrawHisto(signal, background, var, c, h, N);
    Grafica->Run();
	return 0;
}
