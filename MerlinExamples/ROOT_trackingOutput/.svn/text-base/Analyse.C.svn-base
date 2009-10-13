// just a few pedagogical examples
// dk 15.1.2009
//
void Analyse(){

	TFile *merlin = TFile::Open("merlin.root");

	TCanvas* can = new TCanvas("can","Merlin",600,750);

	can->Divide(2,2);

        // sigma_y and gey along the linac
        // adiabatic damping and the constancy of gey
	TH1F* frame = new TH1F("frame","",0,0,10000);
	frame->SetStats(0);
	frame->SetMinimum(0);
	frame->SetMaximum(23);
	can->cd(1);
	frame->GetXaxis()->SetTitle("z[m]");
	frame->GetXaxis()->SetNdivisions(505);
	frame->GetYaxis()->SetTitle("#gamma#epsilon_{y}[nm] | #sigma_{y}[#mum]");
	frame->Draw();
	mytree_0->SetLineColor(kBlue);
	mytree_0->Draw("gey*10^9:z","","lsame"); 

	mytree_0->SetLineColor(kRed);
	mytree_0->Draw("s_y*10^6:z","","lsame"); 
	TLatex* l=new TLatex(5000,21.5,"#gamma#epsilon_{y}");
	l->DrawClone();
	l->SetY(4);
	l->SetTitle("#sigma_{y}");
	l->DrawClone();


	// sigma_y and beta_y
	//
	TH1F* frame2 = new TH1F("frame2","",0,1000,2000);
	frame2->SetStats(0);
	frame2->SetMinimum(0);
	frame2->SetMaximum(6);
	can->cd(2);
	frame2->GetXaxis()->SetTitle("z[m]");
	frame2->GetXaxis()->SetNdivisions(505);
	frame2->GetYaxis()->SetTitle("#sigma_{y}[#mum] | #beta_{y}/100[m]");
	frame2->DrawClone();
	mytree_0->SetLineColor(kRed);
	mytree_0->Draw("s_y*10^6:z","","lsame"); 
	mytree_0->SetLineColor(6);
	mytree_0->Draw("by/100:z","","lsame"); 
	l->SetX(1600);
	l->SetY(5.2);
	l->SetTitle("#sigma_{y}");
	l->DrawClone();
	l->SetY(1.6);
	l->SetTitle("#beta_{y}");
	l->DrawClone();

	// correlation y yp
	// 
	can->cd(3);
	mytree_0_INITIAL->SetMarkerColor(kRed);
	mytree_0_INITIAL->Draw("y:yp");
	mytree_0_FINAL->SetMarkerColor(kGreen);
	mytree_0_FINAL->Draw("y:yp","","same");
	TText* t = new TText(2e-7,15e-6,"initial bunch");
	t->SetTextColor(kRed);
	t->DrawClone();
	t->SetTextColor(kGreen);
	t->SetY(12e-6);
	t->SetTitle("final bunch");
	t->DrawClone();

	// correlation ct dp
	// S-shape due to wakefields in final bunch
	can->cd(4);
	mytree_0_FINAL->Draw("dp:ct","","profile");
}

