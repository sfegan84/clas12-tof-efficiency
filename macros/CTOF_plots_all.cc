 {

  //TLines for overlay
  l_vzLo = new TLine(-9.0, 0.0, -9.0, 5000.0);
  l_vzHi = new TLine(2.0, 0.0, 2.0, 5000.0);
  l_vzLo->SetLineColor(2);
  l_vzHi->SetLineColor(2);

  c1 = new TCanvas("c1","c1",150,10,600,500);

  h_z_vertex_CD_2->GetXaxis()->SetTitle("z-vertex (cm)");
  h_z_vertex_CD_2->Draw();
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

  h_radia_residual_CD_2->GetXaxis()->SetTitle("Distance (cm)");
  h_radia_residual_CD_2->Draw();
  l_dist->Draw();


  c1d = new TCanvas("c1d","c1d",150,10,600,500);
  l_dist2 = new TLine(30.0, 0.0, 30.0, 1000000.0);
  l_dist2->SetLineColor(2);

  h_radia_CTOF_CND_2->GetXaxis()->SetTitle("Distance (cm)");
  h_radia_CTOF_CND_2->Draw();
  l_dist2->Draw();

  c1e = new TCanvas("c1e","c1e",150,10,600,500);
  l_dist3 = new TLine(35.0, 0.0, 35.0, 1000000.0);
  l_dist3->SetLineColor(2);

  h_path_CTOF_CND_2->GetXaxis()->SetTitle("Distance (cm)");
  h_path_CTOF_CND_2->Draw();
  l_dist3->Draw();

  c1f = new TCanvas("c1f","c1f",150,10,600,500);
  h_beta_mom_cut_2->Draw("colz");

  c1g = new TCanvas("c1g","c1g",150,10,600,500);
  TH1F* h_beta_proj = (TH1F*)h_beta_mom_cut_2->ProjectionY("beta_proj",0,250);
  h_beta_proj->SetTitle("Track Beta");
  h_beta_proj->Draw();

  c1h = new TCanvas("c1h","c1h",150,10,600,500);
  l_piLo = new TLine(-0.1, 0.0, -0.1, 25000.0);
  l_piHi = new TLine(0.2, 0.0, 0.2, 25000.0);
  l_pLo = new TLine(0.702, 0.0, 0.702, 50000.0);
  l_pHi = new TLine(1.077, 0.0, 1.077, 50000.0);
  l_piLo->SetLineColor(2);
  l_piHi->SetLineColor(2);
  l_pLo->SetLineColor(2);
  l_pHi->SetLineColor(2);
  pimmass->Draw();
  pimmass2->Draw("same");
  l_piLo->Draw();
  l_piHi->Draw();

  c1ha = new TCanvas("c1ha","c1ha",150,10,600,500);
  pipmass->Draw();
  pipmass2->Draw("same");
  l_piLo->Draw();
  l_piHi->Draw();

  c1hb = new TCanvas("c1hb","c1hb",150,10,600,500);
  protmass->Draw();
  protmass2->Draw("same");
  l_pLo->Draw();
  l_pHi->Draw();

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
  TH2F *h_TracksCDPaddle_Charge_neg_NDF2 = (TH2F*)_file0->Get("h_TracksCDPaddle_Charge_-1_NDF_2_topology_0");
  h_TracksCDPaddle_Charge_neg_NDF2->Draw("colz");

  c2->cd(2);
  TH2F *h_TrajCDPaddle_Charge_neg_NDF2 = (TH2F*)_file0->Get("h_TrajCDPaddle_Charge_-1_NDF_2_topology_0");
  h_TrajCDPaddle_Charge_neg_NDF2->Draw("colz");


  c2_1 = new TCanvas("c2_1","c2_1",150,10,1200,500);
  c2_1->Divide(2,1);

  c2_1->cd(1);
  TH2F *h_TracksCDPaddle_Charge_pos_NDF2 = (TH2F*)_file0->Get("h_TracksCDPaddle_Charge_1_NDF_2_topology_0");
  h_TracksCDPaddle_Charge_pos_NDF2->Draw("colz");

  c2_1->cd(2);
  TH2F *h_TrajCDPaddle_Charge_pos_NDF2 = (TH2F*)_file0->Get("h_TrajCDPaddle_Charge_1_NDF_2_topology_0");
  h_TrajCDPaddle_Charge_pos_NDF2->Draw("colz");

  c2_2 = new TCanvas("c2_2","c2_2",150,10,1200,500);
  c2_2->Divide(2,1);

  c2_2->cd(1);
  TH2F *h_TracksCDPaddle_Charge_neg_NDF2_1 = (TH2F*)_file0->Get("h_TracksCDPaddle_Charge_-1_NDF_2_topology_1");
  h_TracksCDPaddle_Charge_neg_NDF2_1->Draw("colz");

  c2_2->cd(2);
  TH2F *h_TrajCDPaddle_Charge_neg_NDF2_1 = (TH2F*)_file0->Get("h_TrajCDPaddle_Charge_-1_NDF_2_topology_1");
  h_TrajCDPaddle_Charge_neg_NDF2_1->Draw("colz");


  c2_3 = new TCanvas("c2_3","c2_3",150,10,1200,500);
  c2_3->Divide(2,1);

  c2_3->cd(1);
  TH2F *h_TracksCDPaddle_Charge_pos_NDF2_1 = (TH2F*)_file0->Get("h_TracksCDPaddle_Charge_1_NDF_2_topology_1");
  h_TracksCDPaddle_Charge_pos_NDF2_1->Draw("colz");

  c2_3->cd(2);
  TH2F *h_TrajCDPaddle_Charge_pos_NDF2_1 = (TH2F*)_file0->Get("h_TrajCDPaddle_Charge_1_NDF_2_topology_1");
  h_TrajCDPaddle_Charge_pos_NDF2_1->Draw("colz");

  c2_4 = new TCanvas("c2_4","c2_4",150,10,1200,500);
  c2_4->Divide(2,1);

  c2_4->cd(1);
  TH2F *h_TracksCDPaddle_Charge_neg_NDF2_2 = (TH2F*)_file0->Get("h_TracksCDPaddle_Charge_-1_NDF_2_topology_2");
  h_TracksCDPaddle_Charge_neg_NDF2_2->Draw("colz");

  c2_4->cd(2);
  TH2F *h_TrajCDPaddle_Charge_neg_NDF2_2 = (TH2F*)_file0->Get("h_TrajCDPaddle_Charge_-1_NDF_2_topology_2");
  h_TrajCDPaddle_Charge_neg_NDF2_2->Draw("colz");


  c2_5 = new TCanvas("c2_5","c2_5",150,10,1200,500);
  c2_5->Divide(2,1);

  c2_5->cd(1);
  TH2F *h_TracksCDPaddle_Charge_pos_NDF2_2 = (TH2F*)_file0->Get("h_TracksCDPaddle_Charge_1_NDF_2_topology_2");
  h_TracksCDPaddle_Charge_pos_NDF2_2->Draw("colz");

  c2_5->cd(2);
  TH2F *h_TrajCDPaddle_Charge_pos_NDF2_2 = (TH2F*)_file0->Get("h_TrajCDPaddle_Charge_1_NDF_2_topology_2");
  h_TrajCDPaddle_Charge_pos_NDF2_2->Draw("colz");

  c2_6 = new TCanvas("c2_6","c2_6",150,10,1200,500);
  c2_6->Divide(2,1);

  c2_6->cd(1);
  TH2F *h_TracksCDPaddle_Charge_neg_NDF2_3 = (TH2F*)_file0->Get("h_TracksCDPaddle_Charge_-1_NDF_2_topology_3");
  h_TracksCDPaddle_Charge_neg_NDF2_3->Draw("colz");

  c2_6->cd(2);
  TH2F *h_TrajCDPaddle_Charge_neg_NDF2_3 = (TH2F*)_file0->Get("h_TrajCDPaddle_Charge_-1_NDF_2_topology_3");
  h_TrajCDPaddle_Charge_neg_NDF2_3->Draw("colz");


  c2_7 = new TCanvas("c2_7","c2_7",150,10,1200,500);
  c2_7->Divide(2,1);

  c2_7->cd(1);
  TH2F *h_TracksCDPaddle_Charge_pos_NDF2_3 = (TH2F*)_file0->Get("h_TracksCDPaddle_Charge_1_NDF_2_topology_3");
  h_TracksCDPaddle_Charge_pos_NDF2_3->Draw("colz");

  c2_7->cd(2);
  TH2F *h_TrajCDPaddle_Charge_pos_NDF2_3 = (TH2F*)_file0->Get("h_TrajCDPaddle_Charge_1_NDF_2_topology_3");
  h_TrajCDPaddle_Charge_pos_NDF2_3->Draw("colz");



  c2a = new TCanvas("c2a","c2a",150,10,600,500);

  //Efficiency_0_NDF_2->SetMinimum(0.5);
  //Efficiency_0_NDF_2->Draw("colz");

  TH2F *h_Eff_0_NDF_2 = (TH2F*)h_TracksCDPaddle_Charge_neg_NDF2->Clone("Eff_0_NDF_2");
  h_Eff_0_NDF_2->Divide(h_TrajCDPaddle_Charge_neg_NDF2);
  h_Eff_0_NDF_2->SetTitle("Efficiency Neg NDF 2");
  h_Eff_0_NDF_2->SetMinimum(0.5);
  h_Eff_0_NDF_2->Draw("colz");

  c2b = new TCanvas("c2b","c2b",150,10,600,500);

  //Efficiency_1_NDF_2->SetMinimum(0.5);
  //Efficiency_1_NDF_2->Draw("colz");

  TH2F *h_Eff_1_NDF_2 = (TH2F*)h_TracksCDPaddle_Charge_pos_NDF2->Clone("Eff_1_NDF_2");
  h_Eff_1_NDF_2->Divide(h_TrajCDPaddle_Charge_pos_NDF2);
  h_Eff_1_NDF_2->SetTitle("Efficiency Pos NDF 2");
  h_Eff_1_NDF_2->SetMinimum(0.5);
  h_Eff_1_NDF_2->Draw("colz");


  c2c = new TCanvas("c2c","c2c",150,10,600,500);

  //Efficiency_0_NDF_2->SetMinimum(0.5);
  //Efficiency_0_NDF_2->Draw("colz");

  TH2F *h_Eff_0_NDF_2_1 = (TH2F*)h_TracksCDPaddle_Charge_neg_NDF2_1->Clone("Eff_0_NDF_2_1");
  h_Eff_0_NDF_2_1->Divide(h_TrajCDPaddle_Charge_neg_NDF2_1);
  h_Eff_0_NDF_2_1->SetTitle("Efficiency Neg NDF 2 topology 1");
  h_Eff_0_NDF_2_1->SetMinimum(0.5);
  h_Eff_0_NDF_2_1->Draw("colz");

  c2d = new TCanvas("c2d","c2d",150,10,600,500);

  //Efficiency_1_NDF_2->SetMinimum(0.5);
  //Efficiency_1_NDF_2->Draw("colz");

  TH2F *h_Eff_1_NDF_2_1 = (TH2F*)h_TracksCDPaddle_Charge_pos_NDF2_1->Clone("Eff_1_NDF_2_1");
  h_Eff_1_NDF_2_1->Divide(h_TrajCDPaddle_Charge_pos_NDF2);
  h_Eff_1_NDF_2_1->SetTitle("Efficiency Pos NDF 2 topology 1");
  h_Eff_1_NDF_2_1->SetMinimum(0.5);
  h_Eff_1_NDF_2_1->Draw("colz");


  c2e = new TCanvas("c2e","c2e",150,10,600,500);

  //Efficiency_0_NDF_2->SetMinimum(0.5);
  //Efficiency_0_NDF_2->Draw("colz");

  TH2F *h_Eff_0_NDF_2_2 = (TH2F*)h_TracksCDPaddle_Charge_neg_NDF2_2->Clone("Eff_0_NDF_2_2");
  h_Eff_0_NDF_2_2->Divide(h_TrajCDPaddle_Charge_neg_NDF2_2);
  h_Eff_0_NDF_2_2->SetTitle("Efficiency Neg NDF 2 topology 2");
  h_Eff_0_NDF_2_2->SetMinimum(0.5);
  h_Eff_0_NDF_2_2->Draw("colz");

  c2f = new TCanvas("c2f","c2f",150,10,600,500);

  //Efficiency_1_NDF_2->SetMinimum(0.5);
  //Efficiency_1_NDF_2->Draw("colz");

  TH2F *h_Eff_1_NDF_2_2 = (TH2F*)h_TracksCDPaddle_Charge_pos_NDF2_2->Clone("Eff_1_NDF_2_2");
  h_Eff_1_NDF_2_2->Divide(h_TrajCDPaddle_Charge_pos_NDF2_2);
  h_Eff_1_NDF_2_2->SetTitle("Efficiency Pos NDF 2 topology 2");
  h_Eff_1_NDF_2_2->SetMinimum(0.5);
  h_Eff_1_NDF_2_2->Draw("colz");


  c2g = new TCanvas("c2g","c2g",150,10,600,500);

  //Efficiency_0_NDF_2->SetMinimum(0.5);
  //Efficiency_0_NDF_2->Draw("colz");

  TH2F *h_Eff_0_NDF_2_3 = (TH2F*)h_TracksCDPaddle_Charge_neg_NDF2_3->Clone("Eff_0_NDF_2_3");
  h_Eff_0_NDF_2_3->Divide(h_TrajCDPaddle_Charge_neg_NDF2_3);
  h_Eff_0_NDF_2_3->SetTitle("Efficiency Neg NDF 2 topology 3");
  h_Eff_0_NDF_2_3->SetMinimum(0.5);
  h_Eff_0_NDF_2_3->Draw("colz");

  c2h = new TCanvas("c2h","c2h",150,10,600,500);

  //Efficiency_1_NDF_2->SetMinimum(0.5);
  //Efficiency_1_NDF_2->Draw("colz");

  TH2F *h_Eff_1_NDF_2_3 = (TH2F*)h_TracksCDPaddle_Charge_pos_NDF2_3->Clone("Eff_1_NDF_2_3");
  h_Eff_1_NDF_2_3->Divide(h_TrajCDPaddle_Charge_pos_NDF2_3);
  h_Eff_1_NDF_2_3->SetTitle("Efficiency Pos NDF 2 topology 3");
  h_Eff_1_NDF_2_3->SetMinimum(0.5);
  h_Eff_1_NDF_2_3->Draw("colz");




  TH1F *Eff_0_proj_0[48];
  TH1F *Eff_0_proj_1[48];
  TH1F *Eff_0_proj_2[48];
  TH1F *Eff_0_proj_3[48];


  for(int ii=0;ii<48;ii++){
    Eff_0_proj_0[ii] = (TH1F*)h_Eff_0_NDF_2->ProjectionX((Form("eff_0_%d_0",ii+1)),49-ii,49-ii);
    Eff_0_proj_1[ii] = (TH1F*)h_Eff_0_NDF_2_1->ProjectionX(Form("eff_0_%d_1",ii+1),49-ii,49-ii);
    Eff_0_proj_2[ii] = (TH1F*)h_Eff_0_NDF_2_2->ProjectionX(Form("eff_0_%d_2",ii+1),49-ii,49-ii);
    Eff_0_proj_3[ii] = (TH1F*)h_Eff_0_NDF_2_3->ProjectionX(Form("eff_0_%d_3",ii+1),49-ii,49-ii);
    ////Eff_1_proj[ii] = (TH1F*)Efficiency_1_NDF_2->ProjectionX(Form("eff_1_%d",ii+1),52-ii,52-ii);

    Eff_0_proj_0[ii]->GetXaxis()->SetRangeUser(-40,40);
    Eff_0_proj_1[ii]->GetXaxis()->SetRangeUser(-40,40);
    Eff_0_proj_2[ii]->GetXaxis()->SetRangeUser(-40,40);
    Eff_0_proj_3[ii]->GetXaxis()->SetRangeUser(-40,40);
    //Eff_1_proj[ii]->GetXaxis()->SetRangeUser(-50,50);
  }

  c2_effProj = new TCanvas("c2_effProj","c2_effProj",150,10,1200,1000);
  c2_effProj->Divide(4,4);

  for(int i=0;i<16;i++){
    c2_effProj->cd(i+1);
    
    Eff_0_proj_0[i]->Draw("P*");
    Eff_0_proj_1[i]->SetMarkerColor(2);
    Eff_0_proj_1[i]->Draw("P* SAME");
    Eff_0_proj_2[i]->SetMarkerColor(3);
    Eff_0_proj_2[i]->Draw("P* SAME");
    Eff_0_proj_3[i]->SetMarkerColor(4);
    Eff_0_proj_3[i]->Draw("P* SAME");
  }


  c2_effProj_1 = new TCanvas("c2_effProj_1","c2_effProj_1",150,10,1200,1000);
  c2_effProj_1->Divide(4,4);

  for(int i=0;i<16;i++){
    c2_effProj_1->cd(i+1);
    
    Eff_0_proj_0[i+16]->Draw("P*");
    Eff_0_proj_1[i+16]->SetMarkerColor(2);
    Eff_0_proj_1[i+16]->Draw("P* SAME");
    Eff_0_proj_2[i+16]->SetMarkerColor(3);
    Eff_0_proj_2[i+16]->Draw("P* SAME");
    Eff_0_proj_3[i+16]->SetMarkerColor(4);
    Eff_0_proj_3[i+16]->Draw("P* SAME");
  }

  c2_effProj_2 = new TCanvas("c2_effProj_2","c2_effProj_2",150,10,1200,1000);
  c2_effProj_2->Divide(4,4);

  for(int i=0;i<16;i++){
    c2_effProj_2->cd(i+1);
    
    Eff_0_proj_0[i+32]->Draw("P*");
    Eff_0_proj_1[i+32]->SetMarkerColor(2);
    Eff_0_proj_1[i+32]->Draw("P* SAME");
    Eff_0_proj_2[i+32]->SetMarkerColor(3);
    Eff_0_proj_2[i+32]->Draw("P* SAME");
    Eff_0_proj_3[i+32]->SetMarkerColor(4);
    Eff_0_proj_3[i+32]->Draw("P* SAME");
  }

//   c3 = new TCanvas("c3","c3",150,10,1200,500);

//   c3->Divide(2,1);

//   c3->cd(1);
//   TH2F *h_TracksCDPaddle_Charge_neg_NDF3 = (TH2F*)_file0->Get("h_TracksCDPaddle_Charge_-1_NDF_3_topology_2");
//   h_TracksCDPaddle_Charge_neg_NDF3->Draw("colz");

//   c3->cd(2);
//   TH2F *h_TrajCDPaddle_Charge_neg_NDF3 = (TH2F*)_file0->Get("h_TrajCDPaddle_Charge_-1_NDF_3_topology_2");
//   h_TrajCDPaddle_Charge_neg_NDF3->Draw("colz");

//   c3_1 = new TCanvas("c3_1","c3_1",150,10,1200,500);
//   c3_1->Divide(2,1);

//   c3_1->cd(1);
//   TH2F *h_TracksCDPaddle_Charge_pos_NDF3 = (TH2F*)_file0->Get("h_TracksCDPaddle_Charge_1_NDF_3_topology_2");
//   h_TracksCDPaddle_Charge_pos_NDF3->Draw("colz");

//   c3_1->cd(2);
//   TH2F *h_TrajCDPaddle_Charge_pos_NDF3 = (TH2F*)_file0->Get("h_TrajCDPaddle_Charge_1_NDF_3_topology_2");
//   h_TrajCDPaddle_Charge_pos_NDF3->Draw("colz");


//   c3a = new TCanvas("c3a","c3a",150,10,600,500);
//   TH2F *h_Eff_0_NDF_3 = (TH2F*)h_TracksCDPaddle_Charge_neg_NDF3->Clone("Eff_0_NDF_3");
//   h_Eff_0_NDF_3->Divide(h_TrajCDPaddle_Charge_neg_NDF3);
//   h_Eff_0_NDF_3->SetTitle("Efficiency Neg NDF 3");
//   h_Eff_0_NDF_3->SetMinimum(0.5);
//   h_Eff_0_NDF_3->Draw("colz");

//   c3b = new TCanvas("c3b","c3b",150,10,600,500);
//   TH2F *h_Eff_1_NDF_3 = (TH2F*)h_TracksCDPaddle_Charge_pos_NDF3->Clone("Eff_1_NDF_3");
//   h_Eff_1_NDF_3->Divide(h_TrajCDPaddle_Charge_pos_NDF3);
//   h_Eff_1_NDF_3->SetTitle("Efficiency Pos NDF 3");
//   h_Eff_1_NDF_3->SetMinimum(0.5);
//   h_Eff_1_NDF_3->Draw("colz");

//   c4 = new TCanvas("c4","c4",150,10,1200,500);
//   c4->Divide(2,1);

//   c4->cd(1);
//   TH2F *h_TracksCDPaddle_Charge_neg_NDF4 = (TH2F*)_file0->Get("h_TracksCDPaddle_Charge_-1_NDF_4_topology_2");
//   h_TracksCDPaddle_Charge_neg_NDF4->Draw("colz");

//   c4->cd(2);
//   TH2F *h_TrajCDPaddle_Charge_neg_NDF4 = (TH2F*)_file0->Get("h_TrajCDPaddle_Charge_-1_NDF_4_topology_2");
//   h_TrajCDPaddle_Charge_neg_NDF4->Draw("colz");

//   c4_1 = new TCanvas("c4_1","c4_1",150,10,1200,500);
//   c4_1->Divide(2,1);

//   c4_1->cd(1);
//   TH2F *h_TracksCDPaddle_Charge_pos_NDF4 = (TH2F*)_file0->Get("h_TracksCDPaddle_Charge_1_NDF_4_topology_2");
//   h_TracksCDPaddle_Charge_pos_NDF4->Draw("colz");

//   c4_1->cd(2);
//   TH2F *h_TrajCDPaddle_Charge_pos_NDF4 = (TH2F*)_file0->Get("h_TrajCDPaddle_Charge_1_NDF_4_topology_2");
//   h_TrajCDPaddle_Charge_pos_NDF4->Draw("colz");


//   c4a = new TCanvas("c4a","c4a",150,10,600,500);
//   TH2F *h_Eff_0_NDF_4 = (TH2F*)h_TracksCDPaddle_Charge_neg_NDF4->Clone("Eff_0_NDF_4");
//   h_Eff_0_NDF_4->Divide(h_TrajCDPaddle_Charge_neg_NDF4);
//   h_Eff_0_NDF_4->SetTitle("Efficiency Neg NDF 4");
//   h_Eff_0_NDF_4->SetMinimum(0.5);
//   h_Eff_0_NDF_4->Draw("colz");

//   c4b = new TCanvas("c4b","c4b",150,10,600,500);
//   TH2F *h_Eff_1_NDF_4 = (TH2F*)h_TracksCDPaddle_Charge_pos_NDF4->Clone("Eff_1_NDF_4");
//   h_Eff_1_NDF_4->Divide(h_TrajCDPaddle_Charge_pos_NDF4);
//   h_Eff_1_NDF_4->SetTitle("Efficiency Pos NDF 4");
//   h_Eff_1_NDF_4->SetMinimum(0.5);
//   h_Eff_1_NDF_4->Draw("colz");

//   c5 = new TCanvas("c5","c5",150,10,1200,500);
//   c5->Divide(2,1);

//   c5->cd(1);
//   TH2F *h_TracksCDPaddle_Charge_neg_NDF5 = (TH2F*)_file0->Get("h_TracksCDPaddle_Charge_-1_NDF_5_topology_2");
//   h_TracksCDPaddle_Charge_neg_NDF5->Draw("colz");

//   c5->cd(2);
//   TH2F *h_TrajCDPaddle_Charge_neg_NDF5 = (TH2F*)_file0->Get("h_TrajCDPaddle_Charge_-1_NDF_5_topology_2");
//   h_TrajCDPaddle_Charge_neg_NDF5->Draw("colz");


//   c5_1 = new TCanvas("c5_1","c5_1",150,10,1200,500);
//   c5_1->Divide(2,1);

//   c5_1->cd(1);
//   TH2F *h_TracksCDPaddle_Charge_pos_NDF5 = (TH2F*)_file0->Get("h_TracksCDPaddle_Charge_1_NDF_5_topology_2");
//   h_TracksCDPaddle_Charge_pos_NDF5->Draw("colz");

//   c5_1->cd(2);
//   TH2F *h_TrajCDPaddle_Charge_pos_NDF5 = (TH2F*)_file0->Get("h_TrajCDPaddle_Charge_1_NDF_5_topology_2");
//   h_TrajCDPaddle_Charge_pos_NDF5->Draw("colz");


// //   c5a = new TCanvas("c5a","c5a",150,10,600,500);
// //   Efficiency_0_NDF_5->SetMinimum(0.5);
// //   Efficiency_0_NDF_5->Draw("colz");

// //   c5b = new TCanvas("c5b","c5b",150,10,600,500);
// //   Efficiency_1_NDF_5->SetMinimum(0.5);
// //   Efficiency_1_NDF_5->Draw("colz");

//   c6 = new TCanvas("c6","c6",150,10,1200,500);
//   c6->Divide(2,1);

//   c6->cd(1);
//   TH2F *h_TracksCDPaddle_Charge_neg_NDF6 = (TH2F*)_file0->Get("h_TracksCDPaddle_Charge_-1_NDF_6_topology_2");
//   h_TracksCDPaddle_Charge_neg_NDF6->Draw("colz");

//   c6->cd(2);
//   TH2F *h_TrajCDPaddle_Charge_neg_NDF6 = (TH2F*)_file0->Get("h_TrajCDPaddle_Charge_-1_NDF_6_topology_2");
//   h_TrajCDPaddle_Charge_neg_NDF6->Draw("colz");

//   c6_1 = new TCanvas("c6_1","c6_1",150,10,1200,500);
//   c6_1->Divide(2,1);

//   c6_1->cd(1);
//   TH2F *h_TracksCDPaddle_Charge_pos_NDF6 = (TH2F*)_file0->Get("h_TracksCDPaddle_Charge_1_NDF_6_topology_2");
//   h_TracksCDPaddle_Charge_pos_NDF6->Draw("colz");

//   c6_1->cd(2);
//   TH2F *h_TrajCDPaddle_Charge_pos_NDF6 = (TH2F*)_file0->Get("h_TrajCDPaddle_Charge_1_NDF_6_topology_2");
//   h_TrajCDPaddle_Charge_pos_NDF6->Draw("colz");


// //   c6a = new TCanvas("c6a","c6a",150,10,600,500);
// //   Efficiency_0_NDF_6->SetMinimum(0.5);
// //   Efficiency_0_NDF_6->Draw("colz");

// //   c6b = new TCanvas("c6b","c6b",150,10,600,500);
// //   Efficiency_1_NDF_6->SetMinimum(0.5);
// //   Efficiency_1_NDF_6->Draw("colz");

//   c7 = new TCanvas("c7","c7",150,10,1200,500);
//   c7->Divide(2,1);

//   c7->cd(1);
//   TH2F *h_TracksCDPaddle_Charge_neg_NDF7 = (TH2F*)_file0->Get("h_TracksCDPaddle_Charge_-1_NDF_7_topology_2");
//   h_TracksCDPaddle_Charge_neg_NDF7->Draw("colz");


//   c7->cd(2);
//   TH2F *h_TrajCDPaddle_Charge_neg_NDF7 = (TH2F*)_file0->Get("h_TrajCDPaddle_Charge_-1_NDF_7_topology_2");
//   h_TrajCDPaddle_Charge_neg_NDF7->Draw("colz");

//   c7_1 = new TCanvas("c7_1","c7_1",150,10,1200,500);
//   c7_1->Divide(2,1);

//   c7_1->cd(1);
//   TH2F *h_TracksCDPaddle_Charge_pos_NDF7 = (TH2F*)_file0->Get("h_TracksCDPaddle_Charge_1_NDF_7_topology_2");
//   h_TracksCDPaddle_Charge_pos_NDF7->Draw("colz");

//   c7_1->cd(2);
//   TH2F *h_TrajCDPaddle_Charge_pos_NDF7 = (TH2F*)_file0->Get("h_TrajCDPaddle_Charge_1_NDF_7_topology_2");
//   h_TrajCDPaddle_Charge_pos_NDF7->Draw("colz");


// //   c7a = new TCanvas("c7a","c7a",150,10,600,500);
// //   Efficiency_0_NDF_7->SetMinimum(0.5);
// //   Efficiency_0_NDF_7->Draw("colz");

// //   c7b = new TCanvas("c7b","c7b",150,10,600,500);
// //   Efficiency_1_NDF_7->SetMinimum(0.5);
// //   Efficiency_1_NDF_7->Draw("colz");

//   c8 = new TCanvas("c8","c8",150,10,1200,500);
//   c8->Divide(2,1);

//   c8->cd(1);
//   TH2F *h_TracksCDPaddle_Charge_neg_NDF8 = (TH2F*)_file0->Get("h_TracksCDPaddle_Charge_-1_NDF_8_topology_2");
//   h_TracksCDPaddle_Charge_neg_NDF8->Draw("colz");

//   c8->cd(2);
//   TH2F *h_TrajCDPaddle_Charge_neg_NDF8 = (TH2F*)_file0->Get("h_TrajCDPaddle_Charge_-1_NDF_8_topology_2");
//   h_TrajCDPaddle_Charge_neg_NDF8->Draw("colz");

//   c8_1 = new TCanvas("c8_1","c8_1",150,10,1200,500);
//   c8_1->Divide(2,1);

//   c8_1->cd(1);
//   TH2F *h_TracksCDPaddle_Charge_pos_NDF8 = (TH2F*)_file0->Get("h_TracksCDPaddle_Charge_1_NDF_8_topology_2");
//   h_TracksCDPaddle_Charge_pos_NDF8->Draw("colz");

//   c8_1->cd(2);
//   TH2F *h_TrajCDPaddle_Charge_pos_NDF8 = (TH2F*)_file0->Get("h_TrajCDPaddle_Charge_1_NDF_8_topology_2");
//   h_TrajCDPaddle_Charge_pos_NDF8->Draw("colz");


//   c8a = new TCanvas("c8a","c8a",150,10,600,500);
//   Efficiency_0_NDF_8->SetMinimum(0.5);
//   Efficiency_0_NDF_8->Draw("colz");

//   c8b = new TCanvas("c8b","c8b",150,10,600,500);
//   Efficiency_1_NDF_8->SetMinimum(0.5);
//   Efficiency_1_NDF_8->Draw("colz");


}
