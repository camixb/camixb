#include "TLorentzVector.h"
#include "TMath.h"
#include "TRandom3.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TString.h"
#include "TRandom.h"
#include <iostream>
#include "TRandom1.h"
#include "TStyle.h" 
#include <iomanip>
#include <vector>
#include <utility>
#include "particle.h"
#include "TTree.h"


using namespace std;

int main (){
    vector <Evento> dati;  //il primo è theta e il secondo è phi
    
    double m_pi=0.140;//Gev
    double m_K= 0.500; //Gev
    double m_B=5.279; //Gev
    double p_B=0.300; //Gev
    double p_f;
    TLorentzVector p4_B;
    p4_B.SetPxPyPzE(p_B, 0, 0, sqrt(p_B*p_B+m_B*m_B));
    p_f=sqrt(TMath::Power(m_B, 4)+TMath::Power(m_pi, 4)+TMath::Power(m_K, 4)-2*m_B*m_B*m_pi*m_pi-2*m_B*m_B*m_K*m_K-2*m_K*m_K*m_pi*m_pi)/(2*m_B);
    //p_f è il modulo del momento di pi e k nel centro di massa
    
    int N=10000; //numero eventi
    int i;



    int nDau=2;
    double nmass[nDau], p[nDau], theta[nDau], phi[nDau]; 
    TTree* tree = new TTree("datatree", "tree containing our data");
    tree->Branch("p_B", &p_B,  "impulso di B");
    tree->Branch("nDau", &nDau, "prodotti");
    tree->Branch("nmass[Dau]", nmass, "Massa di pi, k");
    tree->Branch("p[Dau]", p, "impulso di pi, k");
    tree->Branch("theta[Dau]", theta, "Theta di pi, k");
    tree->Branch("phi[Dau]", phi, "phi di pi, k");




    vector<double> m_inv;
    TRandom*  gen = new TRandom();
    gen->SetSeed(0);
    double x,y,z;
    for(i=0;i<N;i++){
        //la direzione delle particelle è random
        
        gen->Sphere(x,y,z,p_f);
        TVector3 v_pi(x,y,z);   //creo un 3-vettore. Phi è l'angolo in xy, theta è l'angolo con z
        TVector3 v_k=-v_pi; 
        dati.push_back(Evento(v_pi.Phi(), v_pi.Theta(), v_k.Phi(), v_k.Theta()));
     
    }
    vector<pair<TLorentzVector, TLorentzVector>> cdm; //il primo è pione, il secondo kaone
    vector<pair<TLorentzVector, TLorentzVector>> lab;
    vector<double> angle;
    TLorentzVector pi;
    TLorentzVector k;
    double p_pi_0;
    double p_K_0;
    vector<double> m_inv_meas;
    double p_pi_meas[4];
    double p_k_meas[4];
    double sigma[4]={0.01,0.03, 0.05, 0.1};
    double m_d[4];

    TH1F h_pi("p_pi_meas", "distribution of p_pi", 200, 0, 10);
    TH1F m_1("m_inv of detector #1", "distribution of invariant mass detector #1", 200, 4,10);
    TH1F m_2("m_inv of detector #2", "distribution of invariant mass detector #2", 200, 4,10);
    TH1F m_3("m_inv of detector #3", "distribution of invariant mass detector #3", 200, 4,10);
    TH1F m_4("m_inv of detector #4", "distribution of invariant mass detector #4", 200, 4,10);



    
    
    



    for(vector<Evento>::iterator it=dati.begin(); it!=dati.end(); it++){
        
        pi.SetPxPyPzE(p_f*sin(it->get_pi().second)*cos(it->get_pi().first), p_f*sin(it->get_pi().second)*sin(it->get_pi().first), p_f*cos(it->get_pi().second), sqrt(m_pi*m_pi+p_f*p_f));
        k.SetPxPyPzE(p_f*sin(it->get_k().second)*cos(it->get_k().first), p_f*sin(it->get_k().second)*sin(it->get_k().first), p_f*cos(it->get_k().second), sqrt(m_K*m_K+p_f*p_f));

        cdm.push_back(make_pair(pi, k));
        pi.Boost(p4_B.BoostVector());
        k.Boost(p4_B.BoostVector());

        lab.push_back(make_pair(pi,k));
        m_inv.push_back((k+pi).M());
        angle.push_back(k.Vect().Angle(pi.Vect()));
        p_pi_0=pi.Vect().Mag();
        p_K_0= k.Vect().Mag();

        theta[0]=it->get_pi().second;
        theta[1]=it->get_k().second;
        phi[0]=it->get_pi().first;
        phi[1]=it->get_k().first;
        nmass[0]=m_pi;
        nmass[1]=m_K;
        p[0]=p_pi_0;
        p[1]=p_K_0;
        tree->Fill();
        
        //p_pi_meas_n=gen->Gaus(p_pi_0, sigma*p_pi_0);
        //p_k_meas_n=gen->Gaus(p_K_0, sigma*p_K_0);
        //pi.SetVectMag(pi.Vect(), p_pi_meas);
        //k.SetVectMag(k.Vect(), p_k_meas);
        //m_inv_meas.push_back((k+pi).M());
    
        for(int i=0; i<4; i++){
            p_pi_meas[i]=gen->Gaus(p_pi_0, sigma[i]*p_pi_0);   //per ogni detector vedo il valore misurato di p per pi e K
            p_k_meas[i]=gen->Gaus(p_K_0, sigma[i]*p_K_0);
            pi.SetVectMag(pi.Vect(), p_pi_meas[i]); 
            k.SetVectMag(k.Vect(), p_k_meas[i]);
            m_d[i]=(k+pi).M();
        }
        m_1.Fill(m_d[0]);
        m_2.Fill(m_d[1]);
        m_3.Fill(m_d[2]);
        m_4.Fill(m_d[3]);



    }

    int nbins = 200;
    double mhi= *max_element(m_inv.begin(), m_inv.end());
    double mlo= *min_element(m_inv.begin(), m_inv.end());
    
    TH1F h_m_inv("m_inv", "distribution of invariant mass", nbins, 4, 10);
    h_m_inv.GetXaxis()->SetTitle("Distribution of invariant mass [GeV]");
    double a_hi= *max_element(angle.begin(), angle.end());
    double a_lo= *min_element(angle.begin(), angle.end());
    
    TH1F h_angles("angles", "distribution of opening angles", nbins, a_lo-0.3, a_hi+0.3);
    h_angles.GetXaxis()->SetTitle("Distribution of opening angles");

    //mhi= *max_element(m_inv_meas.begin(), m_inv_meas.end());
    //mlo= *min_element(m_inv_meas.begin(), m_inv_meas.end());
    
    //TH1F h_m_inv_meas("m_inv_meas", "distribution of measured invariant mass", nbins, mlo, mhi);
    //h_m_inv.GetXaxis()->SetTitle("Distribution of measured invariant mass [GeV]");
    
    
    for(int i=0; i<N; i++){
        h_m_inv.Fill(m_inv[i]);
        h_angles.Fill(angle[i]);
     
        //h_m_inv_meas.Fill(m_inv_meas[i]);
        //cout << m_inv[i]<< endl;

    }
    TString rootfname("./data.root");
    TFile* orootfile = new TFile( rootfname, "RECREATE");
     if( !orootfile->IsOpen() ) {
        std::cout << "problems creating root file: exiting... " << std::endl;
        exit(-1);
    }
    std::cout << "storing output in root file " << rootfname << std::endl;



      // Actually write tree in file on disk
    tree->Write();

  // Print some info about the tree
    tree->Print();
    orootfile->Close();

  // ==== Done storing data and using the generator
    delete orootfile;
  // ==== Done using the TTree
    delete tree;
      
    TCanvas canv("canv", "canvas for plotting", 1280, 1024);
    gStyle->SetOptStat(111111);

    h_m_inv.Draw();
    canv.SaveAs("true-mass.pdf");

    h_angles.Draw();
    canv.SaveAs("opening-angle.pdf");

    //h_m_inv_meas.Draw();
    //canv.SaveAs("measured-mass.pdf");

    h_pi.Draw();
    canv.SaveAs("p_pi.pdf");

    h_m_inv.SetMarkerStyle(kFullCircle);
    h_m_inv.Draw();
    m_2.SetMarkerStyle(kFullTriangleDown);
    m_2.Draw("SAME");
    gPad->BuildLegend();
    canv.SaveAs("mass.pdf");

    canv.Divide(2,2);
    canv.cd(1);
    m_1.SetFillColor(kRed);
    m_1.Draw();
    canv.cd(2);
    m_2.SetFillColor(kBlue);
    m_2.Draw();
    canv.cd(3);
    m_3.SetFillColor(kGreen);
    m_3.Draw();
    canv.cd(4);
    m_4.SetFillColor(kYellow);
    m_4.Draw();
    canv.SaveAs("detector.pdf");



    delete gen;
   
   
    return 0;
}