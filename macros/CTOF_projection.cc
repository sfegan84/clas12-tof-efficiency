{

  c1 = new TCanvas("c1","c1",150,10,600,500);
  Efficiency_0_NDF_2->SetMinimum(0.5);
  Efficiency_0_NDF_2->Draw("colz");


  TH1F *Eff_0_proj[48];
  TH1F *Eff_1_proj[48];

  for(int ii=0;ii<48;ii++){
    Eff_0_proj[ii] = (TH1F*)Efficiency_0_NDF_2->ProjectionX(Form("eff_0_%d",ii+1),52-ii,52-ii);
    Eff_1_proj[ii] = (TH1F*)Efficiency_1_NDF_2->ProjectionX(Form("eff_1_%d",ii+1),52-ii,52-ii);

    Eff_0_proj[ii]->GetXaxis()->SetRangeUser(-50,50);
    Eff_1_proj[ii]->GetXaxis()->SetRangeUser(-50,50);
  }

  c2 = new TCanvas("c2","c2",150,10,1200,1000);
  c2->Divide(4,4);
  for(int i=0;i<16;i++){
    //52 -> 37
    c2->cd(i+1);
    Eff_0_proj[i]->Draw();
  }
  
  c3 = new TCanvas("c3","c3",150,10,1200,1000);
  c3->Divide(4,4);
  for(int i=0;i<16;i++){
    c3->cd(i+1);
    //36 -> 21
    Eff_0_proj[i+16]->Draw();
  }
  
  c4 = new TCanvas("c4","c4",150,10,1200,1000);
  c4->Divide(4,4);

  for(int i=0;i<16;i++){
    //20 -> 5
    c4->cd(i+1);
    Eff_0_proj[i+32]->Draw();
  }
  

  c1a = new TCanvas("c1a","c1a",150,10,600,500);
  Efficiency_1_NDF_2->SetMinimum(0.5);
  Efficiency_1_NDF_2->Draw("colz");


  c2a = new TCanvas("c2a","c2a",150,10,1200,1000);
  c2a->Divide(4,4);
  for(int i=0;i<16;i++){
    //52 -> 37
    c2a->cd(i+1);
    Eff_1_proj[i]->Draw();
  }
  
  c3a = new TCanvas("c3a","c3a",150,10,1200,1000);
  c3a->Divide(4,4);
  for(int i=0;i<16;i++){
    c3a->cd(i+1);
    //36 -> 21
    Eff_1_proj[i+16]->Draw();
  }
  
  c4a = new TCanvas("c4a","c4a",150,10,1200,1000);
  c4a->Divide(4,4);

  for(int i=0;i<16;i++){
    //20 -> 5
    c4a->cd(i+1);
    Eff_1_proj[i+32]->Draw();
  }

}
