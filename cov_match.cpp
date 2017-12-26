#include <string.h>
#include <string>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <vector>
#include <cmath>

using namespace std;

class hv_stats
{
// 1 object for each contig seq, statistics for v and h part
public:
	int vlength; //length in bp for viral part
	int hlength; //length in bp for human part
	double vmean; //mean coverage for viral part
	double hmean; //mean coverage for human part
	double vstd; //std coverage for viral part
	double hstd; //std coverage for human part
};

class cov_ratio 
{
// for storing coverage value and ratio
public:
	double cov;
	double ratio;
};

string space_2dot(string text)
{
    replace(text.begin(), text.end(), ' ', '.');
    replace(text.begin(), text.end(), '_', '.');
    return text;
}

void remove_from_string( string &str, char* charsToRemove ) 
{
   	for ( unsigned int i = 0; i < strlen(charsToRemove); ++i ) 
	{
      		str.erase( remove(str.begin(), str.end(), charsToRemove[i]), str.end() );
   	}
}

int main(int argc, char* argv[])
{
// argv[1] string input file hmm output, unix-sorted
// argv[2] string input file genomecov beg -bg, unix-sorted

	string argv1=argv[1];
	string argv2=argv[2];

	ifstream hmm_out(argv1.c_str());
	ifstream gcov(argv2.c_str());
	ofstream output("bed_cov.out");	


	string line;
	for(;getline(hmm_out,line);)
	{
		hv_stats hv;
		hv.vlength=0;hv.hlength=0;hv.vmean=0;hv.hmean=0;hv.vstd=0;hv.hstd=0;

		string seqID;
		istringstream linestream(line);
		getline(linestream,seqID,' ');
		
		string temp;
		vector<string> frag_class; //classification of each 1000 length fragment
		vector<string>::iterator stringi;

		for(;getline(linestream, temp, ' ');)
		{
			if(temp=="viral") hv.vlength=hv.vlength+1000;
			else hv.hlength=hv.hlength+1000;
			frag_class.push_back(temp);
		}
		

		if(hv.vlength==0||hv.hlength==0) continue;

		//read bed file and match, both are sorted
		vector<cov_ratio> hcovrec;
		vector<cov_ratio> vcovrec;

		for(;getline(gcov,line);)
		{
			istringstream tabs(line);
			string col;
			getline(tabs,col,'\t');

			string form_ID=space_2dot(col);
			if(form_ID!=seqID) continue;

			getline(tabs,col,'\t');
			int start=stoi(col)+1; //col originally 0-based, covert to 1-based
			getline(tabs,col,'\t');
			int end=stoi(col);
			getline(tabs,col,'\t');
			double coverage=stod(col);
			
			//recording bed line information

			int start_frag=(start-1)/1000;
			
			int seqendbool=0;
			int frag_num=start_frag;
			for(;(frag_num+1)*1000<end;frag_num++)
			{
				if(frag_num>=frag_class.size()) {seqendbool=1;goto nex;}
				cov_ratio newcov;
				newcov.ratio=(frag_num+1)*1000-max(start,frag_num*1000+1)+1;
				newcov.cov=coverage;

				if(frag_class[frag_num]=="viral") vcovrec.push_back(newcov);
				else hcovrec.push_back(newcov);
			}
			if(frag_num>=frag_class.size()) {seqendbool=1;goto nex;}
			cov_ratio newcov;
			newcov.cov=coverage;
			newcov.ratio=end-max(start,frag_num*1000+1)+1;

			if(frag_class[frag_num]=="viral") vcovrec.push_back(newcov);
			else hcovrec.push_back(newcov);

			nex:
			if(seqendbool==1&&(int)end/1000>=frag_class.size()) break;
		}

		//calculate stats based on bed information
			
		vector<cov_ratio>::iterator covi;
		for(covi=hcovrec.begin();covi<hcovrec.end();covi++)
		{
			hv.hmean=hv.hmean+(double)(*covi).cov*(*covi).ratio/hv.hlength;
		}
		for(covi=hcovrec.begin();covi<hcovrec.end();covi++)
		{
			hv.hstd=hv.hstd+(double)pow(((*covi).cov-hv.hmean),2.0)*(*covi).ratio/hv.hlength;
		}
		hv.hstd=pow(hv.hstd,0.5);

		for(covi=vcovrec.begin();covi<vcovrec.end();covi++)
		{
			hv.vmean=hv.vmean+(double)(*covi).cov*(*covi).ratio/hv.vlength;
		}

		for(covi=vcovrec.begin();covi<vcovrec.end();covi++)
		{
			hv.vstd=hv.vstd+(double)pow(((*covi).cov-hv.vmean),2.0)*(*covi).ratio/hv.vlength;
		}
		hv.vstd=pow(hv.vstd,0.5);

		output<<seqID<<","<<hv.hlength<<","<<hv.vlength<<","<<hv.hmean<<","<<hv.vmean<<","<<hv.hstd<<","<<hv.vstd<<",";
		for(stringi=frag_class.begin();stringi<frag_class.end();stringi++)
		{
			output<<(*stringi)<<"_";
		}
		output<<endl;
	}

return 0;
}
















