/* 
 * File:   AnaData.cc
 * Author: rafopar
 *
 * Created on August 10, 2023, 12:22 PM
 */

#include <cstdlib>
#include <iostream>

// ===== Hipo headers =====
#include <reader.h>
#include <writer.h>
#include <dictionary.h>
#include <vector>
#include <math.h>
using namespace std;

/*
 * 
 */
int main(int argc, char** argv) {





    if (argc <= 2) {
        std::cout << " *** please provide the inputFile name and outputFile name ..." << std::endl;
        cout << "Exiting ..." << endl;
        exit(0);
    }

    auto inputFile=argv[1];
    auto outputFile=argv[2];

    


    hipo::reader reader;
    reader.open(inputFile);

    hipo::dictionary factory;

    reader.readDictionary(factory);

    factory.show();

    hipo::event event;

    hipo::writer writer;
    writer.addDictionary(factory);
    writer.open(outputFile);

    int evCounter = 0;

    
    hipo::bank bRecPart(factory.getSchema("REC::Particle"));

    hipo::bank bFT(factory.getSchema("REC::ForwardTagger"));

    
    int accepted_events=0;
    try {

        while (reader.next() == true) {
            reader.read(event);

            evCounter = evCounter + 1;

            //if( evCounter > 155000 ){break;}
            if (evCounter % 10000 == 0) {
                cout.flush() << "Processed " << evCounter << " events \r";
            }

            event.getStructure(bRecPart);
            event.getStructure(bFT);

            
	    //count the number of p, K-, pi+
	    int n_prot = 0;
            int n_pip = 0;
            int n_Km = 0;

            int nPart = bRecPart.getRows();
	    for (int i_part = 0; i_part < nPart; i_part++) {
	      int pid = bRecPart.getInt("pid", i_part);
	      if (pid == 2212) {
		n_prot+=1;
	      } else if (pid == 211) {
		n_pip+=1;
	      } else if (pid == -321) {
		n_Km+=1;
	      }
	    }
	    if (n_prot==0 or n_pip==0 or n_Km==0)
	      continue;
	    //cout << "K-pi+p found" <<endl;
	    //now check for electron in FT with p>7
	    int n_e_FT=0;
	    for (int i_ft = 0; i_ft < bFT.getRows(); i_ft++){
	      int index = bFT.getInt("pindex", i_ft);
	      if(bRecPart.getInt("pid",index)!=11)
		continue;
	      float p = sqrt(pow(bRecPart.getFloat("px", index),2)
			     +pow(bRecPart.getFloat("py", index),2)
			     +pow(bRecPart.getFloat("pz", index),2));
	      
	      if (10.6 - p>8.0){
		n_e_FT+=1;
	      }
	    }
	    if (n_e_FT>0){
	      writer.addEvent(event);
	      accepted_events+=1;
	      //cout << "event!"<<endl;
	    }
        }
	cout << "input events=" <<evCounter << ";  accepted events=" << accepted_events << endl;
    } catch (const char msg) {
        cerr << msg << endl;
    }


    
    writer.close();
    //writer.showSummary();
    
    return 0;
}
