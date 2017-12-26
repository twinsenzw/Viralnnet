#include <string.h>
#include <string>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <vector>


using namespace std;

void remove_from_string( string &str, char* charsToRemove ) 
{
   	for ( unsigned int i = 0; i < strlen(charsToRemove); ++i ) 
	{
      		str.erase( remove(str.begin(), str.end(), charsToRemove[i]), str.end() );
   	}
}

void linearize(string ifile, string ofile)
{
	string linr="awk -f linearize.awk < ";
	linr=linr + ifile + " | awk -F '\t' '{ if (length($2) >= 10) print}' >" + ifile + "_t";
	system(linr.c_str());

	//getting rid of nonnt symbols
	string ifile2=ifile+"_t";
	ifstream if2(ifile2.c_str());
	ofstream of(ofile.c_str());
	string ID;
	string seq;
	string line;

	for(;getline(if2,line);)
	{
		istringstream linestream(line);
		getline(linestream,ID,'\t');
		getline(linestream,seq,'\t');
		remove_from_string(seq, "nryswkmbdhvuNRYSWKMBDHV.-U");
		of<<ID<<'\t'<<seq<<endl;
	}
	
}

string space_2dot(string text)
{
    replace(text.begin(), text.end(), ' ', '.');
    replace(text.begin(), text.end(), '_', '.');
    return text;
}


void fragmentize(string seq, int win_size, int step, string originalID, ofstream &output, int &count) //count=total # of seqs in the new fasta
{
	
	//remove_from_string(seq, "nryswkmbdhvuNRYSWKMBDHV.-U");
	remove_from_string(originalID, ",");
	for(int i=0; i+win_size<seq.length(); i+=step)
	{
		string window=seq.substr(i,win_size);
		output<<">"<<originalID<<"_subseq_"<<i+1<<"_"<<i+win_size<<endl;
		output<<window<<endl;	
		++count;
	}
}

int main(int argc, char* argv[]) 
// argv[1] - input fasta file
// argv[2] - window size
// argv[3] - step size
{
	vector<string> file_name;
	vector<string>::iterator stringi;
	
	int window_size=stoi(argv[2]);
	int step_size=stoi(argv[3]);

	file_name.push_back(argv[1]);
	for(stringi=file_name.begin();stringi<file_name.end();stringi++) //open a fasta in a list
	{

		/* string	frag_fasta_name	fragmented fasta name 
		   ofstream	frag_fasta	fragmented fasta file
		   int		count		# of seqs in frag_fasta
		   seqsID in frag_fasta in form: >original.ID_subseq_startbp_endbp
		   all ' ' and '_' in originalID replaced by '.'
		*/
		
		int count=0;
		string file_name_linr;
		file_name_linr=*stringi+"_lin.fa";
		linearize(*stringi,file_name_linr);
		
		ifstream contigs_lin(file_name_linr.c_str());
		string line;
		
		string frag_fasta_name=*stringi+"_frag.fa";
		ofstream frag_fasta(frag_fasta_name.c_str());
		for(;getline(contigs_lin,line);) //read sequence from fasta, fragmentize, append to frag fasta file
		{
			istringstream linestream(line);
	
			string originalID, fasta_seq;
			getline(linestream,originalID,'\t');
			originalID=originalID.erase(0,1);
			originalID=space_2dot(originalID);
			getline(linestream,fasta_seq,'\t');

			fragmentize(fasta_seq, window_size, step_size, originalID, frag_fasta, count);
		}
		cout<<count<<endl;
	}



			


return 0;
}
