#include <string.h>
#include <string>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <vector>


using namespace std;

int main(int argc, char* argv[])
{
// argv[1] - string raw.csv filename
// argv[2] - string r initial-classification output filename
// argv[3] - double typeI error
// argv[4] - double typeII error
	

// The construct: 
	string arg1=argv[1];
	string arg2=argv[2];
	string arg3=argv[3];
	string arg4=argv[4];

	double typeI=stod(arg3);
	double typeII=stod(arg4);

	ifstream raw_csv(arg1.c_str());
	ifstream r_output(arg2.c_str());

	/* first write all HMM_states*/

	ofstream HMM_states("HMM_states.csv");

	string line;
	getline(raw_csv,line);
	
	string pres_seqID="";
	vector<string> obs_state;
	vector<string>::iterator stringi;

	int first=1;
	int virus=0;
	int human=0;
	string host_seqID;
	for(;getline(raw_csv,line);)
	{
		istringstream linestream(line);
		host_seqID.clear();
		getline(linestream,host_seqID,'_');
		
		
		if(host_seqID!=pres_seqID) //end of fragments from present host seq
		{
			if(!first) 
			{
				if(virus!=0&&human!=0&&obs_state.size()>1)
				{
					HMM_states<<*(obs_state.begin());

					// HMM prior transition probability: h-h, h-v, v-h, v-v
					//HMM_states<<","<<(double)(human-1.0)/human<<","<<(double)1.0/human<<","<<(double)1.0/virus<<","<<(double)(virus-1.0)/virus;
					HMM_states<<","<<0.999<<","<<0.001<<","<<0.1<<","<<0.9;
					
					// HMM prior emission probability: h-h, h-v, v-h, v-v
					HMM_states<<","<<1-typeI<<","<<typeI<<","<<typeII<<","<<1-typeII;

					HMM_states<<","<<obs_state.size();

					for(stringi=obs_state.begin()+1;stringi<obs_state.end();stringi++)
					{
						HMM_states<<","<<*stringi;
					}
					HMM_states<<endl;
				}
			}
			else first=0;
			
			virus=0;human=0;
			obs_state.clear();
			obs_state.push_back(host_seqID);
			pres_seqID.clear();pres_seqID=host_seqID;
			
		}

		getline(r_output,line);
		obs_state.push_back(line);
		if(line=="viral") virus++;
		else human++;
	
	}

	//-output last contig-	
	HMM_states<<*(obs_state.begin());

	// HMM prior transition probability: h-h, h-v, v-h, v-v
	//HMM_states<<","<<(double)(human-1.0)/human<<","<<(double)1.0/human<<","<<(double)1.0/virus<<","<<(double)(virus-1.0)/virus;
	HMM_states<<","<<0.999<<","<<0.001<<","<<0.1<<","<<0.9;
					
	// HMM prior emission probability: h-h, h-v, v-h, v-v
	HMM_states<<","<<1-typeI<<","<<typeI<<","<<typeII<<","<<1-typeII;

	HMM_states<<","<<obs_state.size();

	for(stringi=obs_state.begin()+1;stringi<obs_state.end();stringi++)
	{
		HMM_states<<","<<*stringi;
	}
	HMM_states<<endl;

	
	
return 0;
			
}
