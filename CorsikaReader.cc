#include "/net/minus/home/shchegolev/CORSIKA/corsika-75600/include/crsRead/MCorsikaReader.h"
#include "/net/minus/home/shchegolev/CORSIKA/corsika-75600/include/crs/TSubBlock.h"
#include "/net/minus/home/shchegolev/CORSIKA/corsika-75600/include/crs/MRunHeader.h"
#include "/net/minus/home/shchegolev/CORSIKA/corsika-75600/include/crs/MEventHeader.h"
#include "/net/minus/home/shchegolev/CORSIKA/corsika-75600/include/crs/MEventEnd.h"
#include "/net/minus/home/shchegolev/CORSIKA/corsika-75600/include/crs/MParticleBlock.h"
#include "/net/minus/home/shchegolev/CORSIKA/corsika-75600/include/crs/MLongitudinalBlock.h"
#include "/net/minus/home/shchegolev/CORSIKA/corsika-75600/include/crs/MParticle.h"

#include <mutex>
#include <thread>
#include <cmath>
#include "generator.h"
#include <iostream>
#include <sstream>
#include <map>
#include <vector>
#include <stdlib.h>
#include <time.h>
using namespace std;

// to hold data of one observation level
struct ObsLevel {
  double x;
  double y;
  double x2;
  double y2;
  double w;
};
struct particle {
    unsigned short int id;
    int p;
    int x;
    int y;
    int t;
};

particle part;
vector < vector <particle> > AllEvent;
vector <particle> OneEvent;
vector <int> ShowerEnergy;
vector <float> ShowerTheta;
vector <float> ShowerPhi;
int ShowerSimCounter=0;
int N;
ofstream fout;
ifstream fin;
string fname;
int OneShowerEnergy;
float OneShowerTheta;
float OneShowerPhi;
std::vector<double> all;
std::mutex g_lock;

void loop(int number, unsigned int seed){
  Generator par2(15, 2, 16, 1.46, 6.0, seed);
  for(unsigned int k=0; k < number; ++k){
      for(int i=0;i<AllEvent.size();i++){
          ++ShowerSimCounter;
          par2.randomcp();
          for(int j=0;j<AllEvent[i].size();j++){
	    par2.set(&AllEvent[i][j].id,&AllEvent[i][j].p,&AllEvent[i][j].x,&AllEvent[i][j].y, &AllEvent[i][j].t);
            par2.process();
          }
	  par2.output(&ShowerSimCounter,&ShowerEnergy[i], &ShowerTheta[i], &ShowerPhi[i]);
      }
  }
    g_lock.lock();
    par2.print_all();
    g_lock.unlock();
}

std::vector<std::thread> thread_pool;

main (int argc, char **argv) 
{  

  cout<<"Input number of throws:";
  cin>>N;

  fout.open("cadrS.dat");
  fout.close();

  fout.open("EASfileOutput.dat");
  fout << "EAS number\tEnergy\tTheta\tPhi\tNumber of particles\n";

  if (argc<2) {
    fin.open("e01.inp");
    getline(fin,fname,'\n');
    fin.close();
  }
  else fname=string(argv[1]);

  crsRead::MCorsikaReader cr(fname, 3);
  
  string inputFile;
  if (fname.rfind('/')==string::npos) {
    inputFile = fname;
  } else {
    inputFile = fname.substr(fname.rfind('/')+1, 
                             fname.size()-fname.rfind('/'));
  }
  
  
  int ShowerCounter = 0;

  crs::MRunHeader Run;
  while (cr.GetRun (Run)) {
    
    crs::MEventHeader Shower;
    while (cr.GetShower(Shower)) {
      

      OneShowerEnergy=int(Shower.GetEnergy());
      OneShowerTheta=Shower.GetTheta();
      OneShowerPhi=Shower.GetPhi();
      
      
      crs::TSubBlock Data;
      while (cr.GetData (Data)) {
        
        switch (Data.GetBlockType ()) {
          
            case crs::TSubBlock::ePARTDATA:
            {
              const crs::MParticleBlock& ParticleData = Data;
              crs::MParticleBlock::ParticleListConstIterator iEntry;
              for (iEntry = ParticleData.FirstParticle();
                   iEntry != ParticleData.LastParticle();
                   ++iEntry) {

		// DUMP
		//iEntry->Dump();
                
                if (iEntry->IsParticle()) {
                  
                  crs::MParticle iPart(*iEntry);
                  
                  part.id    = iPart.GetParticleID();
                  part.p  = int(iPart.GetKinEnergy()*1000);
                  part.x  = int(iPart.GetX());
                  part.y  = int(iPart.GetY());
		  part.t = int(iPart.GetTime()); // time[ns]

                  if ((part.x*part.x+part.y*part.y)<72250000 && part.p>1){
		 OneEvent.push_back(part);
		//out1 << part.id << ' ' << part.p << ' ' << part.x << ' ' << part.y << ' ' << part.t << '\n';
		}	
                }
                
              } // end particle loop
              
              break;
            }
            
            case crs::TSubBlock::eLONG:
              break;
              
            default:
              break;
        } // end data block
        
      } // loop data
      
      crs::MEventEnd ShowerSummary;
      if(cr.GetShowerSummary(ShowerSummary)){
        ++ShowerCounter;
        cout << "---------------------------------\n"
           << " Shower info:\n"
           << "  Shower number = " << ShowerCounter << endl
           <<"Number of particles\t"<< OneEvent.size() <<endl
           << " Primary Energy " << OneShowerEnergy/1000 <<" TeV"<<  endl
           << " EAS Theta " << OneShowerTheta <<  endl
           << " EAS Phi " << OneShowerPhi <<  endl;
        fout << ShowerCounter<<"\t"<<OneShowerEnergy<<"\t"
            << OneShowerTheta<<"\t"<<OneShowerPhi<<"\t"
            <<OneEvent.size()<<endl;
        AllEvent.push_back(OneEvent);
        OneEvent.clear();
        ShowerEnergy.push_back(OneShowerEnergy);
        ShowerTheta.push_back(OneShowerTheta);
        ShowerPhi.push_back(OneShowerPhi);
      }
      else { cout << "---------------------------------\n"
                  << " Shower info:\n"
                  << "  Shower number = " << ShowerCounter+1 << "\n"
                  << "EAS not finished \n";
          OneEvent.clear();}
    } // loop shower

  } // loop runs (usually just 1)
  cout << "---------------------------------\n";
  cout<<"TOTAL NUMBER\t"<<ShowerEnergy.size()<<"\n";

  fout.close();

int number_of_threads = 32;
int step = N/number_of_threads; // std::thread::hardware_concurrency()

for(unsigned int i = 0; i < number_of_threads; ++i){
      if(i==number_of_threads-1) step = N - (number_of_threads-1)*step;
      thread_pool.emplace_back(loop, step, i);

}

for(int i = 0; i < number_of_threads; ++i){
      std::thread& entry = thread_pool[i];
      entry.join();
}

  AllEvent.clear();
  
  return 0;
}


  
