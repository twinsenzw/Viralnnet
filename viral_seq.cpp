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

// argv[1] - string input viral_t-test.out
// argv[2] - string input linearized fast (lin.fa)
// argv[3] - double t-test p cut-off

	string argv1=argv[1];
	string argv2=argv[2];
	string argv3=argv[3];
	string argv2s=argv2+".sorted";
	string argv1c=argv1+".cutoff";

	
	//construct and input cut-off viral_t-test
	string cut_t;
	cut_t="awk '$2 < "+argv3+" ' viral_t-test.out > "+argv1c;
	system(cut_t.c_str());
	ifstream t_test(argv1c.c_str());

	//construct and input sorted linearized fasta
	string sort_lin;
	sort_lin="sort "+argv2+" > "+argv2s;
	system(sort_lin.c_str());
	ifstream lin(argv2s.c_str());

	ofstream vseq("viral_seq.fa");

	string line;
	string line1;
	for(;getline(t_test,line);)
	{
		istringstream linestream(line);
		string seq_ID;
		getline(linestream,seq_ID,' ');
		
		string seq_class;
		getline(linestream,seq_class,' ');
		getline(linestream,seq_class,' ');
		
		string viral_seq; //viral_seq from linearized fasta;
		for(;getline(lin,line1);)
		{
			istringstream linestream1(line1);
			string seq_ID_search;
			getline(linestream1,seq_ID_search,'\t');
			getline(linestream1,viral_seq,'\t');
			seq_ID_search=seq_ID_search.erase(0,1);
			seq_ID_search=space_2dot(seq_ID_search);
			remove_from_string(seq_ID_search, ",");
				
			if(seq_ID_search==seq_ID) break;
		}
		
		istringstream seqstream(seq_class);
		
		string frag;
		int new_seq_bool=1; //start of a new viral sequence?
		int seq_num=0;
		int cnt=0;
		for(;getline(seqstream,frag,'_');cnt++)
		{
			
			if(frag=="viral")
			{
				if(new_seq_bool==1)
				{
					seq_num++;
					vseq<<'>'<<seq_ID<<".seq"<<seq_num<<endl;
					new_seq_bool=0;
				}
				vseq<<viral_seq.substr(cnt*1000,1000);
			}
			else
			{
				if(new_seq_bool==0)
				{
					vseq<<endl;
					new_seq_bool=1;
				}
			}
		}
		if(new_seq_bool==0) vseq<<endl;
	}


return 0;
}











