#include <string>
#include <iostream>
#include <fstream>
#include <cstdio>
#include <ctime>
#include <sstream>
#include <stdlib.h> 
#include <vector>

using namespace std;


const int N = 20, M = 600;
const int L = 11, D = 3;
const string dnaAlphabet = "ACGT";

vector<string> _seqs;
string _motif (L, 'a');
string _editedMotif[N];
int _pos[N];

void generateSeqs(){
for (int i = 0; i < N; ++i)
{
	string temp (M, 'a');
	for (int j = 0; j < M; ++j)
	{
		temp[j] = dnaAlphabet[rand() % dnaAlphabet.size()];
	}
	_seqs.push_back(temp);
}
//for (vector<string>::const_iterator i = _seqs.begin(); i != _seqs.end(); ++i)
  //  std::cout << *i << endl;
}

void generateMotif()
{
	for (int i = 0; i < L; ++i)
	{
		_motif[i] = dnaAlphabet[rand() % dnaAlphabet.size()];
	}
}

string editMotif(string motif)
{
	int i, pos, delta, alpha, beta;
	char x;

	delta = rand()%(D + 1);
	for (i = 0; i < delta; ++i)
	{
		pos = rand()% motif.size();
		motif.erase(pos, 1);
	}

	alpha = rand()%(D - delta + 1);
	for (i = 0; i < alpha; ++i)
	{
		pos = rand()% (motif.size()+1);
		x = dnaAlphabet[rand() % dnaAlphabet.size()];
		motif.insert(pos, 1, x);
	}

	beta = D - delta - alpha;
	for (i = 0; i < beta; ++i)
	{
		pos = rand() % (motif.size());
		x = motif[pos];
		while (x == motif[pos])
		{
			x = dnaAlphabet[rand() % dnaAlphabet.size()];
		}
		motif[pos] = x;
	}
	return motif;
}

void plantedSeqs()
{	
	int endPos;
	
 	for (int i = 0; i < N; ++i)
 	{
		_editedMotif[i] = editMotif(_motif);
		endPos = M - _editedMotif[i].size();
		_pos[i] = rand()%(endPos + 1);
		_seqs[i].replace(_pos[i], _editedMotif[i].size(), _editedMotif[i]);
 	}
}

int main(int argc, char**argv) 
{
	srand ((unsigned int)time(NULL));
	generateSeqs();
	generateMotif();
	plantedSeqs();

	stringstream s;
	s << "planted_l" << L << "_d" << D << ".txt";
	string fileName = s.str();
	ofstream myFile;
    myFile.open(fileName.c_str());
	if (myFile.is_open())
	{
	
		for (int i = 0; i < N; i++)
		{
			myFile << i << " Motif " << _motif << " planted as " << _editedMotif[i] << " at position " << _pos[i] << endl;
			myFile << _seqs[i] << endl;
		}
		myFile.close();
	}
	else 
		cout << "unable to open file";
	return 0;
}
