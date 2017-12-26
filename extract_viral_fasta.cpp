#include <string.h>
#include <string>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <vector>


using namespace std;
class fa_inf  //information for each viral sequence, used to match frag
{
public:
	int start;
	int end;
	string chr;
};

int main(int argc, char* argv[])
{
	string argv1=argv[1];
	ifstream input("hmmout.txt");
	//ifstream fa(argv1.c_str());
	ofstream output("viral_sequences.fa");
	string line;
	string pr_chr=""; //present chromosome
	int pr_cnt=0; //counter for label in the present chromosome
	fa_inf vfrag;
	string fasta;

	for(;getline(input,line);)
	{
		istringstream linestream(line);
		string label;
		int v_extension_bool=0;
		ifstream fa(argv1.c_str());
		for(;linestream >> label;pr_cnt++)
		{
			if(label=="") { pr_cnt--; continue; }
			if(label=="viral")
			{
				if(v_extension_bool==0) //this viral label is the start of a new viral fragment
				{
					vfrag.chr=pr_chr;
					vfrag.start=pr_cnt*1000+1;
					v_extension_bool=1;
					
				}
			}
			else if(label=="human")
			{
				if(v_extension_bool==1) //this human label is met during the extension of a viral fragment, end the extension and finalize the vfrag object
				{
					vfrag.end=pr_cnt*1000;
					v_extension_bool=0;

					//based on the vfrag object, search and output sequence from _frag.fa
					ostringstream startID, endID;
					startID << vfrag.chr << "_subseq_" << vfrag.start << "_" << vfrag.start+1000-1;
					endID << vfrag.chr << "_subseq_" << vfrag.end-1000+1 << "_" << vfrag.end;

					string seq="";
					cout<<startID.str()<<endl;
					cout<<endID.str()<<endl;
					for(;getline(fa,fasta);)
					{
						
						if(fasta.find(startID.str())!=-1)
						{
							cout<<"found"<<endl;
							for(;getline(fa,fasta);)
							{
								if(fasta=="") continue;
								if(startID.str()==endID.str()) //only the first frag is viral
								{
									cout<<"first frag.."<<endl;
									seq.append(fasta);
								        goto end;
								}
								if(fasta[0]=='>')
								{
									
									if(fasta.find(endID.str())!=-1)
									{
										getline(fa,fasta);
										seq.append(fasta);
										goto end;
									}
									
								}
								else
								{
									seq.append(fasta);
								}
							}
						}
					}
					end:
					output<<">"<<vfrag.chr<<"_"<<vfrag.start<<"_"<<vfrag.end<<endl;
					output<<seq<<endl;
				}
			}
			else //other labels: name of a chromosome, update present chr name and present chr position counter
			{
				pr_cnt=-1;
				pr_chr=label;
			}
		}
	}

return 0;

}
