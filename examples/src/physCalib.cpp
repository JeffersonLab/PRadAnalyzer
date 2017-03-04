//============================================================================//
// Automatically search data files between the input begin and end runs       //
// Analyze these data files for physics calibration                           //
//                                                                            //
// Weizhi Xiong                                                               //
// 12/10/2016                                                                 //
//============================================================================//

#include "physCalib.h"
#include "TSystem.h"
#include "TMath.h"
#include "TVector2.h"
#include <algorithm>
#include <cassert>
#include <cmath>
#include <iomanip>

void Helper()
{
    cout << "usage: physCalib <options> <begin_run> <end_run>" << endl
         << "options:" << endl
         << setw(10) << "-i " << "<input_file_dir>" << endl
         << setw(10) << "-o " << "<output_file_dir>" << endl
         << setw(10) << "-h " << "<show_options>" << endl
         << endl;
    exit(1);
}

int main(int argc, char * argv [])
{
    // determine input files
    int run[2], run_arg;
    run[0] = 0;
    run[1] = 999999;
    run_arg = 0;
    string in_dir = "./"; //define your dst data file folder here
    string out_dir = "./";

    for(int i = 1; i < argc; ++i)
    {
        char* ptr = argv[i];
        if(*(ptr++) == '-') {
            switch(*(ptr++))
            {
            case 'o':
                out_dir = argv[++i];
                break;
            case 'i':
                in_dir = argv[++i];
                break;
            case 'h':
            default:
                Helper();
            }
        } else {
            if(run_arg > 1)
                Helper();
            run[run_arg] = atoi(argv[i]);
            run_arg++;
        }
    }

    auto inputFiles = FindInputFiles(in_dir, run[0], run[1]);
    if (inputFiles.empty()) {
        cout << "did not find any input files in dir " << in_dir
             << " between run " << run[0]
             << " and " << run[1]
             << endl;
        return -1;
    }

    hycal_sys = new PRadHyCalSystem("config/hycal.conf");
    gem_sys = new PRadGEMSystem("config/gem.conf");
    epics = new PRadEPICSystem("config/epics_channels.conf");
    coord_sys = new PRadCoordSystem("database/coordinates.dat");
    det_match = new PRadDetMatch("config/det_match.conf");
    dst_parser = new PRadDSTParser();

    dst_parser->SetMode(0);
    hycal = hycal_sys->GetDetector();
    gem1 = gem_sys->GetDetector("PRadGEM1");
    gem2 = gem_sys->GetDetector("PRadGEM2");

    InitInnerModule();
    InitProfileData();

    if(out_dir.back() != '/')
        out_dir += "/";
    string out_file = "cal" + to_string(run[0]) + "_" + to_string(run[1]) + ".root";
    TFile *f = new TFile((out_dir + out_file).c_str(), "RECREATE");

    InitHistogram();  //initialize global histograms
    HyCalZ = 5817.;  //fDetCoor->GetHyCalZ();

    //------------------------------analyze events--------------------------------//
    for (auto &file : inputFiles)
    {
        dst_parser->OpenInput(file.c_str());

        cout << "Open input file: " << file << endl;

        hycal_sys->ChooseRun(file);

        coord_sys->ChooseCoord(PRadInfoCenter::GetRunNumber());

        int count = 0;
        int beam_energy_ch = epics->GetChannel("MBSY2C_energy");
        PRadBenchMark timer;
        while(dst_parser->Read())
        {
            if(dst_parser->EventType() == PRadDSTParser::Type::event) {
                auto event = dst_parser->GetEvent();

                if (!(event.trigger == PHYS_LeadGlassSum || event.trigger == PHYS_TotalSum))
                    continue;

                count++;
                if (count%PROGRESS_COUNT == 0) {
                    cout <<"------[ ev " << count << " ]---"
                         << "---[ " << timer.GetElapsedTimeStr() << " ]---"
                         << "---[ " << timer.GetElapsedTime()/(double)count << " ms/ev ]------"
                         << "\r" << flush;
                }

                current_event = &event;
                //reconstruct HyCal clusters
                hycal_sys->Reconstruct(event);
                //reconstruct GEM clusters
                gem_sys->Reconstruct(event);

                //get reconstructed hits
                auto &hycal_hits = hycal->GetHits();
                auto &gem1_hits = gem1->GetHits();
                auto &gem2_hits = gem2->GetHits();

                // check cluster size
	            clusterN = hycal_hits.size();
	            if (clusterN >= MAXCLUSTER) continue;

                //transform coordinates to beam frame
                coord_sys->Transform(hycal->GetDetID(), hycal_hits.begin(), hycal_hits.end());
                coord_sys->Transform(gem1->GetDetID(), gem1_hits.begin(), gem1_hits.end());
                coord_sys->Transform(gem2->GetDetID(), gem2_hits.begin(), gem2_hits.end());

                // projection, default projection is from target to hycal surface
                coord_sys->Projection(hycal_hits.begin(), hycal_hits.end());
                coord_sys->Projection(gem1_hits.begin(), gem1_hits.end());
                coord_sys->Projection(gem2_hits.begin(), gem2_hits.end());

                // match hits
                auto matched = det_match->Match(hycal_hits, gem1_hits, gem2_hits);

                // fill gem hits if only 1 cluster is detected
                if(gem1_hits.size() == 1) GEM_ep_x_y->Fill(gem1_hits.front().x, gem1_hits.front().y);
                if(gem2_hits.size() == 1) GEM_ep_x_y->Fill(gem2_hits.front().x, gem2_hits.front().y);

                // combine hycalhits and gemhits if matched
                CombinedHit myhits[clusterN];
                for(int i = 0; i < clusterN; ++i)
                {
                    myhits[i] = hycal_hits[i];
                }

                // set the matched gem position
                for(auto hit : matched)
                {
                    myhits[hit.hycal_idx].x_gem = hit.x;
                    myhits[hit.hycal_idx].y_gem = hit.y;
                    // take the matching flag in
                    myhits[hit.hycal_idx].flag = hit.hycal.flag;
                }

	            //get beam energy
	            eventNumber = event.event_number;
	            h_beam_e->Fill(Ebeam);
	            if (clusterN == 1) {
	                EPAnalyzer(myhits);
	            }
	            else if (clusterN == 2) {
	                MollerAnalyzer(myhits);
	                OffsetAnalyzer(myhits);
	            }
	            else if (clusterN > 2) {
                    MultiHitAnalyzer(myhits, clusterN);
                }

            }
            else if(dst_parser->EventType() == PRadDSTParser::Type::epics) {
                const auto &epics_ev = dst_parser->GetEPICSEvent();
	            // save epics into handler, otherwise get epicsvalue won't work
	            epics->AddEvent(epics_ev);
                // only update beam energy when there is an epics event
                Ebeam = epics_ev.values.at(beam_energy_ch);
                // or it can find the value that is closest to current event
                // Ebeam = epics->FindValue(eventNumber,"MBSY2C_energy");
	        }
        }

        dst_parser->CloseInput();
        cout <<"------[ ev " << count << " ]---"
             << "---[ " << timer.GetElapsedTimeStr() << " ]---"
             << "---[ " << timer.GetElapsedTime()/(double)count << " ms/ev ]------"
             << endl;
        cout << "Analyzed " << file << "."
             << endl;
    }

    //----------------------------------------------------------------------------//

    double maxModuleEnergy = 0;
    for (unsigned int i=0; i<innerModule.size(); i++) {
        if (countFill[i] == 0) continue;
        for (int j=0; j<12; j++) {
            for (int k=0; k<12; k++) {
                //moduleEnergy[i][j][k] /= (float)countFill[i];
                if (moduleEnergy[i][j][k] > maxModuleEnergy) maxModuleEnergy = moduleEnergy[i][j][k];
                profile[i]->Fill(k, j, moduleEnergy[i][j][k]);
            }
        }
        profile[i]->Fill(5, 5, -0.01*maxModuleEnergy);
        profile[i]->Fill(6, 6, -0.01*maxModuleEnergy);
        profile[i]->Fill(5, 6, -0.01*maxModuleEnergy);
        profile[i]->Fill(6, 5, -0.01*maxModuleEnergy);
    }

    f->cd();

    h_beam_e->Write();
    ep_x_y->Write();
    GEM_x_y->Write();
    GEM_ep_x_y->Write();
    sym_ee_x_y->Write();
    ep_ratio_all->Write();
    ep_angle_E->Write();
    ee1_x_y->Write();
    ee2_x_y->Write();
    ee_ratio_all->Write();
    ee_angle_E->Write();
    sym_ee_E->Write();
    eloss->Write();
    total_angle_E->Write();

    TDirectory *cur_dir = f->mkdir("Modules");
    cur_dir->cd();
    for(int i=0;i<T_BLOCKS;i++) {
        if(ep_ratio[i])
            ep_ratio[i]->Write();
        if(ep_energy[i])
            ep_energy[i]->Write();
        if(ee_ratio[i])
            ee_ratio[i]->Write();
        if(ee_energy[i])
            ee_energy[i]->Write();
    }

    cur_dir = f->mkdir("Moller Coords");
    cur_dir->cd();
    for (int i=0; i<2; i++) {
        deltaCoor[i]->Write();
        coorIntersect[i]->Write();
        deltaCoor_GEM[i]->Write();
        coorIntersect_GEM[i]->Write();
        sym_ee_r[i]->Write();
    }

    cur_dir = f->mkdir("Profile");
    cur_dir->cd();
    for (unsigned int i=0; i<profile.size(); i++) {
        profile[i]->Write();
    }

    cur_dir = f->mkdir("Inner Modules");
    cur_dir->cd();
    for (int i=0; i<12 ; i++) {
        inner_ratio[i]->Write();
        inner_deltax[i]->Write();
        inner_deltay[i]->Write();
        ratio_x[i]->Write();
        ratio_y[i]->Write();
    }
  f->Close();

  return 0;
}
//____________________________________________________________________________________________
void InitHistogram()
{
    float Ehigh = 0;
    int nBin = 0;
    if (isLowE) {
        Ehigh = 1500.;
        nBin = 150;
    }
    else{
        Ehigh = 2700;
        nBin = 270;
    }

    h_beam_e = new TH1F("h_beam_e","beam E",500, 0, Ehigh);
    ep_x_y = new TH2F("ep_x_y","ep cluster position",500,-600,600,500,-600,600);
    ep_x_y->SetOption("colz");
    GEM_x_y = new TH2F("GEM_x_y","GEM hiy position",500,-600,600,500,-600,600);
    GEM_x_y->SetOption("colz");
    GEM_ep_x_y = new TH2F("GEM_ep_x_y","GEM ep hit position",500,-600,600,500,-600,600);
    GEM_ep_x_y->SetOption("colz");
    ep_ratio_all = new TH1F("ep_ratio_all","ep ratio all",500,0.0,2.0);
    ep_angle_E = new TH2F("ep_angle_E","ep recE vs angle",500,0.0,8.0,500,0,Ehigh);
    ep_angle_E->SetOption("colz");
    ee1_x_y = new TH2F("ee1_x_y","ee1 cluster position",500,-600,600,500,-600,600);
    ee1_x_y->SetOption("colz");
    ee2_x_y = new TH2F("ee2_x_y","ee2 cluster position",500,-600,600,500,-600,600);
    ee2_x_y->SetOption("colz");
    ee_ratio_all = new TH1F("ee_ratio_all","ee ratio all",500,0.0,2.0);
    ee_angle_E = new TH2F("ee_angle_E","ee recE vs angle",500,0.0,8.0,500,0,Ehigh);
    ee_angle_E->SetOption("colz");
    sym_ee_x_y = new TH2F("sym_ee_x_y", "sym_ee_x_y", 500, -600, 600, 500, -600, 600);
    sym_ee_x_y->SetOption("colz");
    sym_ee_E = new TH1F("sym_ee_E", "sym_ee_E", 500, Ehigh/2. - 500, Ehigh/2. + 500);
    eloss = new TH2F("eloss", "eloss", 500, 0, Ehigh, 500, 0, 5);
    eloss->SetOption("colz");
    total_angle_E = new TH2F("total_angle_E", "total_angle_E", 500, 0, 8.0, 500, 0, Ehigh);
    total_angle_E->SetOption("colz");

    for(int i=0;i<T_BLOCKS;i++) {
        PRadHyCalModule *module = hycal->GetModule(i+1);
        if(module) {
            string ep_name1 = "ep_ratio_" + module->GetName();
            string ep_name2 = "ep_energy_" + module->GetName();
            string ee_name1 = "ee_ratio_" + module->GetName();
            string ee_name2 = "ee_energy_" + module->GetName();
            ep_ratio[i] = new TH1F(ep_name1.c_str(),ep_name1.c_str(),200,0,2);
            ep_energy[i] = new TH1F(ep_name2.c_str(),ep_name2.c_str(),nBin,0,Ehigh);
            ee_ratio[i] = new TH1F(ee_name1.c_str(),ee_name1.c_str(),200,0,2);
            ee_energy[i] = new TH1F(ee_name2.c_str(),ee_name2.c_str(),nBin,0,Ehigh);
        } else {
            ep_ratio[i] = nullptr;
            ep_energy[i] = nullptr;
            ee_ratio[i] = nullptr;
            ee_energy[i] = nullptr;
        }
    }

    for (int i=0; i<2; i++) {
        deltaCoor[i] = new TH1F(Form("delta_coor_%d", i+1), Form("delta_coor_%d", i+1), 500, -100, 100 );
        coorIntersect[i] = new TH1F(Form("intersect_%d", i+1), Form("intersect_%d", i+1), 500, -50, 50);
        deltaCoor_GEM[i] = new TH1F(Form("delta_coor_%d_GEM", i+1), Form("delta_coor_%d_GEM", i+1), 500, -20, 20 );
        coorIntersect_GEM[i] = new TH1F(Form("intersect_%d_GEM", i+1), Form("intersect_%d_GEM", i+1), 500, -20, 20);
        sym_ee_r[i] = new TH1F(Form("sym_ee_r_%d", i+1), Form("sym_ee_r_%d", i+1), 500, 50, 350);
    }

    for (unsigned int i=0; i<innerModule.size(); i++) {
        TH2F* thisPlot = new TH2F(Form("profile_%d", innerModule[i]), Form("profile_%d", innerModule[i]), 12, 0, 12, 12, 0, 12);
        thisPlot->SetOption("LEGO1");
        profile.push_back(thisPlot);
    }

    for (int i=0; i<12; i++) {
        inner_ratio[i] = new TH1F(Form("inner_ratio_%d", innerModuleList[i]), Form("inner_ratio_%d", innerModuleList[i]), 500,  0, 100);
        inner_deltax[i] = new TH1F(Form("inner_deltax_%d", innerModuleList[i]), Form("inner_deltax_%d", innerModuleList[i]), 500,  -60, 60);
        inner_deltay[i] = new TH1F(Form("inner_deltay_%d", innerModuleList[i]), Form("inner_deltay_%d", innerModuleList[i]), 500,  -60, 60);
        ratio_x[i] = new TH2F(Form("ratio_x_%d", innerModuleList[i]), Form("ratio_x_%d", innerModuleList[i]), 1000, -80, 80, 1000, 0, 2);
        ratio_x[i]->SetOption("colz");
        ratio_y[i] = new TH2F(Form("ratio_y_%d", innerModuleList[i]), Form("ratio_y_%d", innerModuleList[i]), 1000, -80, 80, 1000, 0, 2);
        ratio_y[i]->SetOption("colz");
    }
}
//____________________________________________________________________________________________
void InitProfileData()
{
    ConfigParser parser;
    if (!parser.OpenFile("./database/prof_pwo.dat")) {
        cout<<"cannot find profile data for PWO"<<endl;
        return;
    }

    while (parser.ParseLine()) {
        int i = parser.TakeFirst().Int();
        int j = parser.TakeFirst().Int();
        float val = parser.TakeFirst().Float();
        pwo_profile[i][j] = val;
        pwo_profile[j][i] = val;
        parser.TakeFirst().Float();
    }
    parser.CloseFile();
}
//____________________________________________________________________________________________
inline void EPAnalyzer(CombinedHit* h)
{
    if (h[0].nblocks <= 3) return; //some channel has over charge

    float ratio = h[0].E / GetExpectedEnergy(kProton, h[0].x, h[0].y);
    ep_ratio_all->Fill(ratio);
    ep_x_y->Fill(h[0].x, h[0].y);
	ep_angle_E->Fill(atan(sqrt(h[0].x*h[0].x +
	                 h[0].y*h[0].y)/HyCalZ)*180/TMath::Pi(), h[0].E);
	total_angle_E->Fill(atan(sqrt(h[0].x*h[0].x +
	                 h[0].y*h[0].y)/HyCalZ)*180/TMath::Pi(), h[0].E);
	if (h[0].cid > 0) {

	    //GEM histogram if within range

        if (fabs(h[0].x) < innerBoundary[0] && fabs(h[0].y) < innerBoundary[1]) {
            InnerModuleAnalyzer(h);
        }
	    else{
	        ep_ratio[h[0].cid - 1]->Fill(ratio);
	        ep_energy[h[0].cid - 1]->Fill(h[0].E);
	    }
	}

	if (isLowE && h[0].E > 1300) ee2_x_y->Fill(h[0].x, h[0].y);

}
//____________________________________________________________________________________________
void InnerModuleAnalyzer(CombinedHit* h)
{
    if (! MatchedGEM(h[0]) || !TEST_BIT(h[0].flag, kPbWO4)) return;

    // we need module energy, so need to refresh all modules' energies by choose current event
    hycal_sys->ChooseEvent(*current_event);

    //we work here in the internal coordinate system of HyCal for now, unit in cm
    //x axis points to beam left, y points up
    float expectE = GetExpectedEnergy(kProton, h[0].x_gem, h[0].y_gem);

	float gemX = h[0].x_gem;
	float gemY = h[0].y_gem;
	PRadHyCalModule *module = hycal->GetModule(gemX, gemY);
    if(module == nullptr) return;
    int id = module->GetID();
	FillInnerHist(id, h[0].E, expectE, kProton);

	int neighborID = IsNeighborToInner(id);
	if (neighborID > 0) {
	    PRadHyCalModule *thisM = hycal->GetModule(neighborID);
	    float thisX = (thisM->GetX() - gemX)/2.077;
	    float thisY = (thisM->GetY() - gemY)/2.075;
	    float expectPotion = GetExpectedEFromProfile(thisX, thisY);
	    //cout<<expectPotion<<" "<<thisX<<" "<<thisY<<" "<<gemX<<" "<<gemY<<" "<<neighborID<<" "<<id<<endl;
	    if (expectPotion > 0.) {

	        PRadHyCalModule* module = hycal->GetModule(neighborID);
	        if (module != nullptr) {
	            int indexSave = -1;
	            for (int i=0; i<12; i++) {
	                if (neighborID == innerModuleList[i]) {
	                    indexSave = i;
	                    break;
	                }
	            }
	            if (indexSave >= 0) {
	                if (h[0].E > 900) {
	                    inner_ratio[indexSave]->Fill(module->GetEnergy() / (expectE*expectPotion));
	                    inner_deltax[indexSave]->Fill(h[0].x - h[0].x_gem);
	                    inner_deltay[indexSave]->Fill(h[0].y - h[0].y_gem);
	                }
	                ratio_x[indexSave]->Fill(h[0].x_gem, h[0].E/expectE);
	                ratio_y[indexSave]->Fill(h[0].y_gem, h[0].E/expectE);
	            }
	            //if (indexSave >= 0) inner_ratio[indexSave]->Fill( sqrt( pow(h[0].x, 2) + pow(h[0].y, 2) ) -  sqrt( pow(h[0].x_gem, 2) + pow(h[0].y_gem, 2) ));
	        }
	    }
	}

    PRadHyCalModule *m = hycal->GetModule(id);
    if(m == nullptr) return;
    float dx = m->GetX() - gemX;
    float dy = m->GetY() - gemY;
	float dist = sqrt(dx*dx + dy*dy);
	if (dist > 0.3) return;


	if (fabs(h[0].E - Ebeam)/Ebeam > 0.1) return;
	bool pass = true;
	int  indexSave = -1;
	for (unsigned int i=0; i<innerModule.size(); i++) {
	    if (id == innerModule[i]) {
	        pass = false;
	        indexSave = i;
	        break;
	    }
	}
	if (pass) return;


	GEM_x_y->Fill(h[0].x_gem, h[0].y_gem);
	//GEM_x_y->Fill(gemX*10., gemY*-10.);

	countFill[indexSave]++;

	for (unsigned int i=0; i<boundModule.size(); i++) {
	    PRadHyCalModule* module = hycal->GetModule(boundModule[i]);
	    if (module == nullptr) continue;
	    if (module->GetEnergy() < 1.4) continue; //1.4 is the threshold for PWO modules

	    int indexj = (boundModule[i]-1000)%34 - 12;
	    int indexi = 11 - ((boundModule[i]-1000)/34 - 11);

	    moduleEnergy[indexSave][indexi][indexj] += module->GetEnergy();
	}
}
//____________________________________________________________________________________________
inline void FillInnerHist(int& id, float& E, float& expectE, ParticleType type)
{
    if (type == kProton) {
        ep_ratio[id - 1]->Fill(E / expectE);
	    ep_energy[id - 1]->Fill(E);
    }else{
        ee_ratio[id - 1]->Fill(E / expectE);
	    ee_energy[id - 1]->Fill(E);
    }
}
//____________________________________________________________________________________________
inline void MollerAnalyzer(CombinedHit* h, int index1, int index2)
{
    int in[2] = {index1, index2};
    float r[2], phi[2], ratio[2], theta[2];
    for (int i=0; i<2; i++) {
        r[i]   = sqrt(h[in[i]].x*h[in[i]].x + h[in[i]].y*h[in[i]].y);
        phi[i] = h[in[i]].x>0 ? atan(h[in[i]].y / h[in[i]].x)*180.0/TMath::Pi() :
                                    atan(h[in[i]].y / h[in[i]].x)*180.0/TMath::Pi() + 180.;
    }
    //not co-plane enough, not considered as Moller, see if it is likely a ep
    if (fabs( fabs(phi[0] - phi[1]) - 180.) > 20.) {
        for (int i=0; i<2; i++) {
            if( ( (TEST_BIT(h[in[i]].flag, kPbGlass) || TEST_BIT(h[in[i]].flag, kTransition) ) && h[in[i]].E > (1.-0.2727)*Ebeam) ||
                ( (TEST_BIT(h[in[i]].flag, kPbWO4) && h[in[i]].E > (1.-0.0909)*Ebeam) ) ) EPAnalyzer(&h[in[i]]);
        }
        return;
    }

    for (int i=0; i<2; i++) {
        //ratio[i] = h[in[i]].E / weight[i];
        ratio[i] = h[in[i]].E / GetExpectedEnergy(kElectron, h[in[i]].x, h[in[i]].y);
        theta[i] = atan(sqrt(h[in[i]].x*h[in[i]].x + h[in[i]].y*h[in[i]].y)/HyCalZ)*180./TMath::Pi();
        //check if this is actually a EP event
        if ( ( (TEST_BIT(h[in[i]].flag, kPbGlass) || TEST_BIT(h[in[i]].flag, kTransition) ) && h[in[i]].E > (1.-0.2727)*Ebeam)
            || (TEST_BIT(h[in[i]].flag, kPbWO4) && theta[i] > 1.5 && h[in[i]].E > (1.-0.0909)*Ebeam) ) {
            //Moller cannot have such high energy in LG and Transition region
            EPAnalyzer(&h[in[i]]);
            break;
        }
        ee_ratio_all->Fill(ratio[i]);
        ee1_x_y->Fill(h[in[i]].x, h[in[i]].y);
        ee_angle_E->Fill(theta[i], h[in[i]].E);
        total_angle_E->Fill(theta[i], h[in[i]].E);
        if (h[in[i]].cid > 0) {

            //GEM histogram if within range
            if (fabs(h[in[i]].x) < innerBoundary[0] && fabs(h[in[i]].y) < innerBoundary[1]) {
                if (MatchedGEM(h[in[i]]) && TEST_BIT(h[in[i]].flag, kPbWO4)) {
                    float expectE = GetExpectedEnergy(kElectron, h[in[i]].x_gem, h[in[i]].y_gem);
                    float gemX = h[in[i]].x_gem;
	                float gemY = h[in[i]].y_gem;
	                PRadHyCalModule *module = hycal->GetModule(gemX, gemY);
                    if (module) {
                        int id = module->GetID();
                        FillInnerHist(id, h[in[i]].E, expectE, kElectron);
                    }
                }
            }else{
                ee_ratio[h[in[i]].cid-1]->Fill(ratio[i]);
                ee_energy[h[in[i]].cid-1]->Fill(h[in[i]].E);
            }
        }
    }

    //if are symmetric Moller, each has about half of the beam energy, HyCal resolution assumed to be 0.03
    if (! ( fabs( h[index1].E - Ebeam/2.) < 0.03*Ebeam && fabs(h[index2].E - Ebeam/2.) < 0.03*Ebeam) ) return;
    if (clusterN != 2) return;
    if (theta[0] < 0.7 || theta[1] < 0.7) return;
    if (h[index1].x < 0) {
        sym_ee_r[0]->Fill(r[0]); sym_ee_r[1]->Fill(r[1]);
    }else{
        sym_ee_r[0]->Fill(r[1]); sym_ee_r[1]->Fill(r[0]);
    }
    sym_ee_x_y->Fill(h[index1].x, h[index1].y);
    sym_ee_x_y->Fill(h[index2].x, h[index2].y);
    sym_ee_E->Fill(h[index1].E);
    sym_ee_E->Fill(h[index2].E);
}
//____________________________________________________________________________________________
inline void MultiHitAnalyzer(CombinedHit* h, int& n)
{
    if (n >= MAXCLUSTER) return;
    std::sort(h, h+n, SortForE);

    float weight[MAXCLUSTER], phi[MAXCLUSTER];

    for (int i=0; i<n; i++) {
        phi[i] = h[i].x>0 ? atan(h[i].y / h[i].x)*180.0/TMath::Pi() :
                                atan(h[i].y / h[i].x)*180.0/TMath::Pi() + 180.;
    }
    //find the best match for the highest energy cluster
    for (int i=1; i<clusterN ; i++) {
        weight[i] = 0.7 * fabs( (Ebeam - h[0].E - h[i].E) )/30. + 0.3 * fabs( fabs(phi[0] - phi[i]) - 180.) / 10.;
    }

    int idSave = -1;
    float bestWeight = 1e9;
    for (int i=1; i<clusterN; i++) {
        if (bestWeight > weight[i]) {
            bestWeight = weight[i];
            idSave = i;
        }
    }
    MollerAnalyzer(h, 0, idSave);
}
//_________________________________________________________________________________________________
inline void OffsetAnalyzer(CombinedHit* h)
{
    assert(clusterN == 2); //we want only two cluster event to do this

    if (!MatchedGEM(h[0]) || !MatchedGEM(h[1])) return;

    //we apply a set of rather tight cut on Moller candidates that could participate in this calculation
    if ((h[0].E + h[1].E) < (1. - 0.15)*Ebeam || (h[0].E + h[1].E) > (1. + 0.15)*Ebeam) return;
    if (!TEST_BIT(h[0].flag, kPbWO4) || !TEST_BIT(h[1].flag, kPbWO4)) return; //PWO has better position resolution
    double phi[2] = { atan2(h[0].y, h[0].x), atan2(h[1].y, h[1].x) };
    double theta[2] = { atan( sqrt( h[0].x*h[0].x + h[0].y*h[0].y)/HyCalZ )*180./TMath::Pi(),
                        atan( sqrt( h[1].x*h[1].x + h[1].y*h[1].y)/HyCalZ )*180./TMath::Pi() };

    if ( fabs( TVector2::Phi_mpi_pi(phi[1] - phi[0] - TMath::Pi()) ) >  10.*TMath::Pi() / 180.) return;
    if (theta[0] < 0.7 || theta[1] < 0.7) return;

    float xinter, yinter;
    GetIntersect(h[0].x, h[1].x, h[0].y, h[1].y, xinter, yinter);
    coorIntersect[0]->Fill(xinter);
    coorIntersect[1]->Fill(yinter);
    GetIntersect(h[0].x_gem, h[1].x_gem, h[0].y_gem, h[1].y_gem, xinter, yinter);
    coorIntersect_GEM[0]->Fill(xinter);
    coorIntersect_GEM[1]->Fill(yinter);

    assert(pairHyCal.size() < 3);
    if (pairHyCal.size() < 2) {
       pairHyCal.push_back(PairCoor(h[0].x, h[0].y, h[1].x, h[1].y));
       pairGEM.push_back(PairCoor(h[0].x_gem, h[0].y_gem, h[1].x_gem, h[1].y_gem));
    }else{
       float xsect, ysect, xsect_GEM, ysect_GEM;
       LinesIntersect(xsect, ysect, xsect_GEM, ysect_GEM);
       deltaCoor[0]->Fill(xsect);
       deltaCoor[1]->Fill(ysect);
       deltaCoor_GEM[0]->Fill(xsect_GEM);
       deltaCoor_GEM[1]->Fill(ysect_GEM);
       pairHyCal.clear();
       pairGEM.clear();
    }
}
//_________________________________________________________________________________________________
inline void LinesIntersect(float &xsect, float &ysect, float &xsect_GEM, float & ysect_GEM)
{
    float xa[2];
    float ya[2];
    float xb[2];
    float yb[2];

    for (unsigned int i=0; i<pairHyCal.size(); i++) {
        xa[i] = pairHyCal.at(i).x1;
        ya[i] = pairHyCal.at(i).y1;
        xb[i] = pairHyCal.at(i).x2;
        yb[i] = pairHyCal.at(i).y2;
    }

    float m[2];
    float b[2];

    for (int i=0; i<2; i++) {
      m[i] = (yb[i] - ya[i]) / (xb[i] - xa[i]);
      b[i] = m[i]*xa[i] - ya[i];
    }

    xsect = (b[1] - b[0])/(m[0] - m[1]);
    ysect = (m[1]*b[0] - m[0]*b[1])/(m[1] - m[0]);
    //             GEM              //
    for (unsigned int i=0; i<pairGEM.size(); i++) {
        xa[i] = pairGEM.at(i).x1;
        ya[i] = pairGEM.at(i).y1;
        xb[i] = pairGEM.at(i).x2;
        yb[i] = pairGEM.at(i).y2;
    }

    for (int i=0; i<2; i++) {
      m[i] = (yb[i] - ya[i]) / (xb[i] - xa[i]);
      b[i] = m[i]*xa[i] - ya[i];
    }

    xsect_GEM = -1.*(b[1] - b[0])/(m[0] - m[1]);
    ysect_GEM = -1.*(m[1]*b[0] - m[0]*b[1])/(m[1] - m[0]);
}
//____________________________________________________________________________________________
vector<string> FindInputFiles(const string &in_dir, int start, int end)
{
    vector<string> res;
    void* dirp = gSystem->OpenDirectory(in_dir.c_str());
    if (!dirp) {
        cerr << "Error: Can't open input file directory "<< in_dir << endl;
        return res;
    }
    cout << "Scanning for input files in " << in_dir
         << " within run " << start << " ~ " << end << "..." << flush;

    const char* dir_item;
    while( (dir_item = gSystem->GetDirEntry(dirp)) ) {
        string file_name = dir_item;
        if ((file_name.find("prad") != string::npos) &&
            (file_name.substr(file_name.size() - 3) == "dst")) {
            string thisName = in_dir;
            if(in_dir.back() != '/')
                thisName.append("/");
            thisName.append(file_name);

            int runNumber = GetRunNumber(thisName);
            if (runNumber <= end && runNumber >= start) {
                res.push_back(thisName);
            }
        }
    }
    cout << endl << "Found " << res.size() << " input files " << endl;
    std::sort(res.begin(), res.end(), SortFile);
    gSystem->FreeDirectory(dirp);
    return res;

}
//______________________________________________________________________________________________
float GetExpectedEnergy(ParticleType type, float& x, float& y)
{
    float theta = atan(sqrt(x*x + y*y)/HyCalZ);
    float expectE = Ebeam*targetMass[type] / ( Ebeam*(1.-cos(theta)) + targetMass[type] );
    float eLoss = GetElossIonElectron(theta, expectE);
    eloss->Fill(expectE, eLoss);
    return expectE - eLoss;
}
//______________________________________________________________________________________________
float GetExpectedEFromProfile(float& dx, float& dy)
{
    float ax = 100*fabs(dx);
    float ay = 100*fabs(dy);
    int i = (int)ax;
    int j = (int)ay;
    if (i<0 || j<0 || i>499 || j>499) return 0.;
    float wx = ax-i;
    float wy = ay-j;
    return pwo_profile[i][j]*wx*wy+pwo_profile[i+1][j]*(1-wx)*wy+pwo_profile[i][j+1]*wx*(1-wy)+pwo_profile[i+1][j+1]*(1-wx)*(1-wy);
}
//______________________________________________________________________________________________
void GetIntersect(float& x1, float& x2, float& y1, float &y2, float &x, float &y)
{
    float m = (y2 - y1) / (x2 - x1);
    y = -1.*m*x1 + y1;
    x = -1.*y1/m + x1;
}
//______________________________________________________________________________________________
int GetRunNumber(string path)
{
    string name = ConfigParser::decompose_path(path).name;
    return ConfigParser::find_integer(name);
}
//______________________________________________________________________________________________
bool SortFile(string name1, string name2) {
    return (GetRunNumber(name1) < GetRunNumber(name2)) ;
}
//______________________________________________________________________________________________
bool SortForE(CombinedHit const & a, CombinedHit const & b)
{
    return a.E > b.E;
}
//______________________________________________________________________________________________
float GetElossIonElectron(float &theta, float& E)
{
    // Calculates energy loss dE/dx in MeV/mm due to ionization for relativistic electrons/positrons.
    //
    // For formula used, see:
    // Track fitting with energy loss
    // Stampfer, Regler and Fruehwirth
    // Computer Physics Communication 79 (1994), 157-164
    //
    // ZoverA:     atomic number / atomic mass of passed material
    // density:    density of material in g/mm^3
    // I : mean excitation energy in MeV

    //only the Al thin window now, need to add GEM, GEM frame, and the cover between HyCal and GEM to be exact

    float ZoverA[3]   = { 13./27., 10.6/21.8 , 0.49919};
    float density[3]  = {2.699, 0.1117433, 1.205e-3};                      // g/cm^2
    float I[3]        = {166*1.e-6, 106.6e-6, 85.7e-6};                       // MeV
    float de       = 5.0989 * 1.e-25;                 // 4*pi*re*me*c^2 in MeV * mm^2 with re = classical electron radius
    float avogadro = TMath::Na();                     // Avogadro constant in 1/mol
    float me       = 0.5109989181;   // electron mass in MeV/c^2
    float gamma    = E / me;                          // Relativistic gamma-factor.

    // Formula is slightly different for electrons and positrons.
    float gammaFac = 3.;
    float corr     = 1.95;

    float eDep = 0;
    float length[3] = { 0.2, 1.5, 50 };
    for (int i=0; i<3; i++) {
      length[i] /= cos(theta);
      float dedx = 0.5 * de * avogadro * density[i] * ZoverA[i] * (2 * TMath::Log(2*me/I[i]) + gammaFac * TMath::Log(gamma) - corr);
      eDep += dedx*length[i];
    }

    return eDep;//convert from MeV/mm to GeV/m
}
//___________________________________________________________________________________________________
void InitInnerModule()
{
    int startCol = 14;
    int startRow = 13;
    int count = 0;
    for (int i=0; i<8; i++) {
        for (int j=0; j<8; j++) {
            innerModule.push_back(1000 + (startRow+i)*34 + startCol+j);
            countFill.push_back(0);
            count++;
        }
    }
    assert(count == 64);

    startCol = 12;
    startRow = 11;
    count = 0;
    for (int i=0; i<12; i++) {
        for (int j=0; j<12; j++) {
            boundModule.push_back(1000 + (startRow+i)*34 + startCol+j);
            count++;
        }
    }
    assert(count == 144);
}
//______________________________________________________________________________________________________
inline int IsNeighborToInner(int& id)
{
    if (id == 1492) return 1526;
    else if (id == 1493) return 1527;
    else if (id == 1494) return 1528;
    else if (id == 1495) return 1529;
    else if (id == 1530) return 1529;
    else if (id == 1564) return 1563;
    else if (id == 1598) return 1597;
    else if (id == 1632) return 1631;
    else if (id == 1665) return 1631;
    else if (id == 1664) return 1630;
    else if (id == 1663) return 1629;
    else if (id == 1662) return 1628;
    else if (id == 1627) return 1628;
    else if (id == 1593) return 1594;
    else if (id == 1559) return 1560;
    else if (id == 1525) return 1526;
    else return -1;
    /*if (id == 1492) return 1458;
    else if (id == 1493) return 1459;
    else if (id == 1494) return 1460;
    else if (id == 1495) return 1461;
    else if (id == 1530) return 1531;
    else if (id == 1564) return 1565;
    else if (id == 1598) return 1599;
    else if (id == 1632) return 1633;
    else if (id == 1665) return 1699;
    else if (id == 1664) return 1698;
    else if (id == 1663) return 1697;
    else if (id == 1662) return 1696;
    else if (id == 1627) return 1626;
    else if (id == 1593) return 1592;
    else if (id == 1559) return 1558;
    else if (id == 1525) return 1524;
    else return -1;*/

}

bool MatchedGEM(const CombinedHit &hit)
{
    return TEST_BIT(hit.flag, kGEM1Match) || TEST_BIT(hit.flag, kGEM2Match);
}

