{

  //TLines for overlay
  l_vzLo = new TLine(-9.0, 0.0, -9.0, 5000.0);
  l_vzHi = new TLine(2.0, 0.0, 2.0, 5000.0);
  l_vzLo->SetLineColor(2);
  l_vzHi->SetLineColor(2);

  c1 = new TCanvas("c1","c1",150,10,600,500);

  h_z_vertex_CD->GetXaxis()->SetTitle("z-vertex (cm)");
  h_z_vertex_CD->Draw();
  l_vzLo->Draw();
  l_vzHi->Draw();
  

  c1a = new TCanvas("c1a","c1a",150,10,1200,500);
  c1a->Divide(2,1);

  c1a->cd(1);
  h_TracksNDF_CD_0->GetXaxis()->SetRangeUser(0,20);
  h_TracksNDF_CD_0->Draw();

  c1a->cd(2);
  h_TracksNDF_CD_1->GetXaxis()->SetRangeUser(0,20);
  h_TracksNDF_CD_1->Draw();

  c1c = new TCanvas("c1c","c1c",150,10,600,500);
  l_dist = new TLine(6.0, 0.0, 6.0, 1000000.0);
  l_dist->SetLineColor(2);

  h_radia_residual_CD->GetXaxis()->SetTitle("Distance (cm)");
  h_radia_residual_CD->Draw();
  l_dist->Draw();


  c1d = new TCanvas("c1d","c1d",150,10,600,500);
  l_dist2 = new TLine(30.0, 0.0, 30.0, 1000000.0);
  l_dist2->SetLineColor(2);

  h_radia_CTOF_CND->GetXaxis()->SetTitle("Distance (cm)");
  h_radia_CTOF_CND->Draw();
  l_dist2->Draw();

  c1e = new TCanvas("c1e","c1e",150,10,600,500);
  l_dist3 = new TLine(35.0, 0.0, 35.0, 1000000.0);
  l_dist3->SetLineColor(2);

  h_path_CTOF_CND->GetXaxis()->SetTitle("Distance (cm)");
  h_path_CTOF_CND->Draw();
  l_dist3->Draw();

  c1f = new TCanvas("c1f","c1f",150,10,600,500);
  h_beta_mom_cut->Draw("colz");

  c1g = new TCanvas("c1g","c1g",150,10,600,500);
  TH1F* h_beta_proj = (TH1F*)h_beta_mom_cut->ProjectionY("beta_proj",0,250);
  h_beta_proj->SetTitle("Track Beta");
  h_beta_proj->Draw();

  c1h = new TCanvas("c1h","c1h",150,10,600,500);
  //l_piLo = new TLine(-0.1, 0.0, -0.1, 25000.0);
  //l_piHi = new TLine(0.2, 0.0, 0.2, 25000.0);
  l_piLo = new TLine(0.702, 0.0, 0.702, 50000.0);
  l_piHi = new TLine(1.077, 0.0, 1.077, 50000.0);
  l_piLo->SetLineColor(2);
  l_piHi->SetLineColor(2);
  //pimmass->Draw();
  //pimmass2->Draw("same");
  //pipmmass->Draw();
  //pipmmass2->Draw("same");
  prmmass->Draw();
  prmmass2->Draw("same");
  l_piLo->Draw();
  l_piHi->Draw();

  c1i = new TCanvas("c1i","c1i",150,10,600,500);
  DeltaMomentum->Draw();

  c1j = new TCanvas("c1j","c1j",150,10,600,500);
  DeltaTheta->Draw();

  c1k = new TCanvas("c1k","c1k",150,10,600,500);
  DeltaPhi->Draw();

  c1l = new TCanvas("c1l","c1l",150,10,1200,1000);
  c1l->Divide(2,2);
  c1l->cd(1);
  el_thetaPhi->Draw("colz");
  c1l->cd(2);
  prot_thetaPhi->Draw("colz");
  c1l->cd(3);
  pipl_thetaPhi->Draw("colz");
  c1l->cd(4);
  pimi_thetaPhi->Draw("colz");

  c2 = new TCanvas("c2","c2",150,10,1200,500);
  c2->Divide(2,1);

  c2->cd(1);
  TH2F *h_TracksCD_Charge_neg_NDF2 = (TH2F*)_file0->Get("h_TracksCD_Charge_-1_NDF_2");
  h_TracksCD_Charge_neg_NDF2->Draw("colz");

  c2->cd(2);
  TH2F *h_TrajCD_Charge_neg_NDF2 = (TH2F*)_file0->Get("h_TrajCD_Charge_-1_NDF_2");
  h_TrajCD_Charge_neg_NDF2->Draw("colz");


  c2_1 = new TCanvas("c2_1","c2_1",150,10,1200,500);
  c2_1->Divide(2,1);

  c2_1->cd(1);
  TH2F *h_TracksCD_Charge_pos_NDF2 = (TH2F*)_file0->Get("h_TracksCD_Charge_1_NDF_2");
  h_TracksCD_Charge_pos_NDF2->Draw("colz");

  c2_1->cd(2);
  TH2F *h_TrajCD_Charge_pos_NDF2 = (TH2F*)_file0->Get("h_TrajCD_Charge_1_NDF_2");
  h_TrajCD_Charge_pos_NDF2->Draw("colz");


  c2a = new TCanvas("c2a","c2a",150,10,600,500);

  Efficiency_0_NDF_2->SetMinimum(0.5);
  Efficiency_0_NDF_2->Draw("colz");

  c2b = new TCanvas("c2b","c2b",150,10,600,500);

  Efficiency_1_NDF_2->SetMinimum(0.5);
  Efficiency_1_NDF_2->Draw("colz");


  c3 = new TCanvas("c3","c3",150,10,1200,500);

  c3->Divide(2,1);

  c3->cd(1);
  TH2F *h_TracksCD_Charge_neg_NDF3 = (TH2F*)_file0->Get("h_TracksCD_Charge_-1_NDF_3");
  h_TracksCD_Charge_neg_NDF3->Draw("colz");

  c3->cd(2);
  TH2F *h_TrajCD_Charge_neg_NDF3 = (TH2F*)_file0->Get("h_TrajCD_Charge_-1_NDF_3");
  h_TrajCD_Charge_neg_NDF3->Draw("colz");

  c3_1 = new TCanvas("c3_1","c3_1",150,10,1200,500);
  c3_1->Divide(2,1);

  c3_1->cd(1);
  TH2F *h_TracksCD_Charge_pos_NDF3 = (TH2F*)_file0->Get("h_TracksCD_Charge_1_NDF_3");
  h_TracksCD_Charge_pos_NDF3->Draw("colz");

  c3_1->cd(2);
  TH2F *h_TrajCD_Charge_pos_NDF3 = (TH2F*)_file0->Get("h_TrajCD_Charge_1_NDF_3");
  h_TrajCD_Charge_pos_NDF3->Draw("colz");


  c3a = new TCanvas("c3a","c3a",150,10,600,500);
  Efficiency_0_NDF_3->SetMinimum(0.5);
  Efficiency_0_NDF_3->Draw("colz");

  c3b = new TCanvas("c3b","c3b",150,10,600,500);
  Efficiency_1_NDF_3->SetMinimum(0.5);
  Efficiency_1_NDF_3->Draw("colz");

  c4 = new TCanvas("c4","c4",150,10,1200,500);
  c4->Divide(2,1);

  c4->cd(1);
  TH2F *h_TracksCD_Charge_neg_NDF4 = (TH2F*)_file0->Get("h_TracksCD_Charge_-1_NDF_4");
  h_TracksCD_Charge_neg_NDF4->Draw("colz");

  c4->cd(2);
  TH2F *h_TrajCD_Charge_neg_NDF4 = (TH2F*)_file0->Get("h_TrajCD_Charge_-1_NDF_4");
  h_TrajCD_Charge_neg_NDF4->Draw("colz");

  c4_1 = new TCanvas("c4_1","c4_1",150,10,1200,500);
  c4_1->Divide(2,1);

  c4_1->cd(1);
  TH2F *h_TracksCD_Charge_pos_NDF4 = (TH2F*)_file0->Get("h_TracksCD_Charge_1_NDF_4");
  h_TracksCD_Charge_pos_NDF4->Draw("colz");

  c4_1->cd(2);
  TH2F *h_TrajCD_Charge_pos_NDF4 = (TH2F*)_file0->Get("h_TrajCD_Charge_1_NDF_4");
  h_TrajCD_Charge_pos_NDF4->Draw("colz");


  c4a = new TCanvas("c4a","c4a",150,10,600,500);
  Efficiency_0_NDF_4->SetMinimum(0.5);
  Efficiency_0_NDF_4->Draw("colz");

  c4b = new TCanvas("c4b","c4b",150,10,600,500);
  Efficiency_1_NDF_4->SetMinimum(0.5);
  Efficiency_1_NDF_4->Draw("colz");

  c5 = new TCanvas("c5","c5",150,10,1200,500);
  c5->Divide(2,1);

  c5->cd(1);
  TH2F *h_TracksCD_Charge_neg_NDF5 = (TH2F*)_file0->Get("h_TracksCD_Charge_-1_NDF_5");
  h_TracksCD_Charge_neg_NDF5->Draw("colz");

  c5->cd(2);
  TH2F *h_TrajCD_Charge_neg_NDF5 = (TH2F*)_file0->Get("h_TrajCD_Charge_-1_NDF_5");
  h_TrajCD_Charge_neg_NDF5->Draw("colz");


  c5_1 = new TCanvas("c5_1","c5_1",150,10,1200,500);
  c5_1->Divide(2,1);

  c5_1->cd(1);
  TH2F *h_TracksCD_Charge_pos_NDF5 = (TH2F*)_file0->Get("h_TracksCD_Charge_1_NDF_5");
  h_TracksCD_Charge_pos_NDF5->Draw("colz");

  c5_1->cd(2);
  TH2F *h_TrajCD_Charge_pos_NDF5 = (TH2F*)_file0->Get("h_TrajCD_Charge_1_NDF_5");
  h_TrajCD_Charge_pos_NDF5->Draw("colz");


  c5a = new TCanvas("c5a","c5a",150,10,600,500);
  Efficiency_0_NDF_5->SetMinimum(0.5);
  Efficiency_0_NDF_5->Draw("colz");

  c5b = new TCanvas("c5b","c5b",150,10,600,500);
  Efficiency_1_NDF_5->SetMinimum(0.5);
  Efficiency_1_NDF_5->Draw("colz");

  c6 = new TCanvas("c6","c6",150,10,1200,500);
  c6->Divide(2,1);

  c6->cd(1);
  TH2F *h_TracksCD_Charge_neg_NDF6 = (TH2F*)_file0->Get("h_TracksCD_Charge_-1_NDF_6");
  h_TracksCD_Charge_neg_NDF6->Draw("colz");

  c6->cd(2);
  TH2F *h_TrajCD_Charge_neg_NDF6 = (TH2F*)_file0->Get("h_TrajCD_Charge_-1_NDF_6");
  h_TrajCD_Charge_neg_NDF6->Draw("colz");

  c6_1 = new TCanvas("c6_1","c6_1",150,10,1200,500);
  c6_1->Divide(2,1);

  c6_1->cd(1);
  TH2F *h_TracksCD_Charge_pos_NDF6 = (TH2F*)_file0->Get("h_TracksCD_Charge_1_NDF_6");
  h_TracksCD_Charge_pos_NDF6->Draw("colz");

  c6_1->cd(2);
  TH2F *h_TrajCD_Charge_pos_NDF6 = (TH2F*)_file0->Get("h_TrajCD_Charge_1_NDF_6");
  h_TrajCD_Charge_pos_NDF6->Draw("colz");


  c6a = new TCanvas("c6a","c6a",150,10,600,500);
  Efficiency_0_NDF_6->SetMinimum(0.5);
  Efficiency_0_NDF_6->Draw("colz");

  c6b = new TCanvas("c6b","c6b",150,10,600,500);
  Efficiency_1_NDF_6->SetMinimum(0.5);
  Efficiency_1_NDF_6->Draw("colz");

  c7 = new TCanvas("c7","c7",150,10,1200,500);
  c7->Divide(2,1);

  c7->cd(1);
  TH2F *h_TracksCD_Charge_neg_NDF7 = (TH2F*)_file0->Get("h_TracksCD_Charge_-1_NDF_7");
  h_TracksCD_Charge_neg_NDF7->Draw("colz");


  c7->cd(2);
  TH2F *h_TrajCD_Charge_neg_NDF7 = (TH2F*)_file0->Get("h_TrajCD_Charge_-1_NDF_7");
  h_TrajCD_Charge_neg_NDF7->Draw("colz");

  c7_1 = new TCanvas("c7_1","c7_1",150,10,1200,500);
  c7_1->Divide(2,1);

  c7_1->cd(1);
  TH2F *h_TracksCD_Charge_pos_NDF7 = (TH2F*)_file0->Get("h_TracksCD_Charge_1_NDF_7");
  h_TracksCD_Charge_pos_NDF7->Draw("colz");

  c7_1->cd(2);
  TH2F *h_TrajCD_Charge_pos_NDF7 = (TH2F*)_file0->Get("h_TrajCD_Charge_1_NDF_7");
  h_TrajCD_Charge_pos_NDF7->Draw("colz");


  c7a = new TCanvas("c7a","c7a",150,10,600,500);
  Efficiency_0_NDF_7->SetMinimum(0.5);
  Efficiency_0_NDF_7->Draw("colz");

  c7b = new TCanvas("c7b","c7b",150,10,600,500);
  Efficiency_1_NDF_7->SetMinimum(0.5);
  Efficiency_1_NDF_7->Draw("colz");

  c8 = new TCanvas("c8","c8",150,10,1200,500);
  c8->Divide(2,1);

  c8->cd(1);
  TH2F *h_TracksCD_Charge_neg_NDF8 = (TH2F*)_file0->Get("h_TracksCD_Charge_-1_NDF_8");
  h_TracksCD_Charge_neg_NDF8->Draw("colz");

  c8->cd(2);
  TH2F *h_TrajCD_Charge_neg_NDF8 = (TH2F*)_file0->Get("h_TrajCD_Charge_-1_NDF_8");
  h_TrajCD_Charge_neg_NDF8->Draw("colz");

  c8_1 = new TCanvas("c8_1","c8_1",150,10,1200,500);
  c8_1->Divide(2,1);

  c8_1->cd(1);
  TH2F *h_TracksCD_Charge_pos_NDF8 = (TH2F*)_file0->Get("h_TracksCD_Charge_1_NDF_8");
  h_TracksCD_Charge_pos_NDF8->Draw("colz");

  c8_1->cd(2);
  TH2F *h_TrajCD_Charge_pos_NDF8 = (TH2F*)_file0->Get("h_TrajCD_Charge_1_NDF_8");
  h_TrajCD_Charge_pos_NDF8->Draw("colz");


  c8a = new TCanvas("c8a","c8a",150,10,600,500);
  Efficiency_0_NDF_8->SetMinimum(0.5);
  Efficiency_0_NDF_8->Draw("colz");

  c8b = new TCanvas("c8b","c8b",150,10,600,500);
  Efficiency_1_NDF_8->SetMinimum(0.5);
  Efficiency_1_NDF_8->Draw("colz");


}
