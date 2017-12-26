#include <string.h>
#include <string>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <vector>

using namespace std;



int main(int argc,char* argv[])
{
// argv[1] - input csv file name
// argv[2] - output csv file name
// argv[3] - label to append to csv

	string input_name=argv[1];
	string output_name=argv[2];

	ifstream feature(input_name.c_str());
	ofstream feature_prep(output_name.c_str());

	string label=argv[3]; //The label to be added to the end column;
	string line;
	getline(feature,line);
	for(;getline(feature,line);)
	{
		vector<double> kmer_rate;
		vector<double>::iterator doublei;
		
		istringstream linestream(line);
		string temp;
		double sum=0;
		getline(linestream,temp,',');
		for(;getline(linestream,temp,',');)
		{
			double count=stod(temp);
			kmer_rate.push_back(count);
			sum=sum+count;
		}
		feature_prep<<label<<","<<(*kmer_rate.begin())/sum;
		for(doublei=kmer_rate.begin()+1;doublei<kmer_rate.end();++doublei)
		{
			feature_prep<<","<<(*doublei)/sum;
		}
		feature_prep<<endl;
	}

return 0;
}
