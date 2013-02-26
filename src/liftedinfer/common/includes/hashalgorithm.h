#ifndef __LVRHASH_ALGORITHM
#define __LVRHASH_ALGORITHM
#include "lvrmln.h"

struct LvrHashAlgorithm
{
	//format of input for hashing grounded atom: <normparentid,term1,term2,....>
	//Note:Only to be called for atoms that are completely grounded (no variables)
	static int DJBHash(vector<int> vals)
	{
		augmentVector(vals);
		unsigned int hash = 5381;
		for(std::size_t i = 0; i < vals.size(); i++)
		{
			hash = ((hash << 5) + hash) + vals[i];
		}
		return hash;
	}

	static unsigned int DJBHashUS(vector<int> vals)
	{
		augmentVector(vals);
		unsigned int hash = 5381;
		for(std::size_t i = 0; i < vals.size(); i++)
		{
			hash = ((hash << 5) + hash) + vals[i];
		}
		return hash;
	}

	static int convertToHash(Atom* atom)
	{
		if(atom->isConstant())
		{
			vector<int> signature;
			signature.push_back(atom->symbol->normParentId);
			for(unsigned int jj=0;jj<atom->terms.size();jj++)
			{
				signature.push_back(atom->terms[jj]->domain[0]);
			}
			return DJBHash(signature);

		}
		vector<vector<int> > signature;
		vector<int> id(1);
		id[0] = atom->symbol->normParentId;
		signature.push_back(id);
		for(unsigned int i=0;i<atom->terms.size();i++)
		{
			vector<int> trm(atom->terms[i]->domain.size());
			for(unsigned int j=0;j<atom->terms[i]->domain.size();j++)
				trm[j] = atom->terms[i]->domain[j];
			signature.push_back(trm);
		}
		return DJBLHash(signature);
	}

	//format of input for hashing atom with variables: <normparentid><[val(if term is constant)][<domval0><domval1>....(if variable)]>....
	//TO be only called for atoms with atleast 1 variable in its terms
	static int DJBLHash(vector<vector<int> > vals)
	{
		vector<int> signature;
		signature.push_back(vals[0].at(0));
		for(unsigned int i=1;i<vals.size();i++)
		{
			for(unsigned int j=0;j<vals[i].size();j++)
				signature.push_back(vals[i].at(j));
		}
		return DJBHash(signature);
		/*
		unsigned int hash = 5381;
		for(std::size_t i = 0; i < signature.size(); i++)
		{
			hash = ((hash << 5) + hash) + signature[i];
		}
		return hash;
		*/
	}
	
	//augment the vector with positional information
	static void augmentVector(vector<int>& input)
	{
		vector<int> output;
		int pos = 1;
		for(unsigned int i=0;i<input.size();i++)
		{
			output.push_back(pos*input[i]);
			output.push_back(input[i]);
			pos++;
		}
		input.clear();
		input=output;
	}
	
};


#endif
