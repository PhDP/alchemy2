#ifndef __LVRATOMHASH_TEMPLATE
#define __LVRATOMHASH_TEMPLATE
#include "hashalgorithm.h"
#include <map>
using namespace std;

template<typename T>
struct LvrAtomHashTemplate
{
private:
	map<int,T> data;
public:
	LvrAtomHashTemplate(){}
	~LvrAtomHashTemplate(){}

	void insertOrReplaceValue(Atom* atom,T value_)
	{

		int key = LvrHashAlgorithm::convertToHash(atom);
		T value = value_;
		typename map<int,T>::iterator it = data.find(key);
		if(it == data.end())
		{
			data.insert(pair<int,T>(key,value));
		}
		else
		{
			it->second = (T) value;
		}
	}
	bool find(Atom* atom)
	{
		int key = LvrHashAlgorithm::convertToHash(atom);
		if(data.find(key)!=data.end())
			return true;
		else
			return false;
	}
	void update(Atom* atom,T value_)
	{
		int key = LvrHashAlgorithm::convertToHash(atom);
		T value = value_;
		typename map<int,T>::iterator it = data.find(key);
		if(it!=data.end())
		{
			it->second = value;
		}
	}

	void update(vector<int> signature,T value_)
	{
		int key = LvrHashAlgorithm::DJBHash(signature);
		T value = value_;
		typename map<int,T>::iterator it = data.find(key);
		if(it!=data.end())
		{
			it->second = value;
		}
	}

	
	void insert(Atom* atom,T value_)
	{
		int key = LvrHashAlgorithm::convertToHash(atom);
		T value = value_;
		typename map<int,T>::iterator it = data.find(key);
		data.insert(pair<int,T>(key,value));
	}

	//representation: int[0]=atom's normparentid, followed by each term's grounded value
	void insert(vector<int> signature,T value_)
	{
		int key = LvrHashAlgorithm::DJBHash(signature);		
		T value = value_;
		typename map<int,T>::iterator it = data.find(key);
		data.insert(pair<int,T>(key,value));
	}

	void incrementValue(vector<int> signature,int value_)
	{
		if(value_==0)
			return;
		int key = LvrHashAlgorithm::DJBHash(signature);
		int value = value_;
		map<int,int>::iterator it = data.find(key);
		if(it!=data.end())
		{
			it->second += value;
		}
	}

	void insertOrincrementValue(Atom* atom,int value_)
	{
		if(value_==0)
			return;
		int key = LvrHashAlgorithm::convertToHash(atom);
		int value = value_;
		map<int,int>::iterator it = data.find(key);
		data.insert(pair<int,int>(key,value));
		if(it==data.end())
		{
			it->second += value;
		}
	}

	void incrementValue(Atom* atom,int value_)
	{
		if(value_==0)
			return;
		int key = LvrHashAlgorithm::convertToHash(atom);
		int value = value_;
		map<int,int>::iterator it = data.find(key);
		if(it!=data.end())
		{
			it->second += value;
		}
	}

	T getValue(Atom* atom)
	{
		//convert to int rep
		int key = LvrHashAlgorithm::convertToHash(atom);
		typename map<int,T>::iterator it = data.find(key);
		if(it!=data.end())
		{
			return it->second;
		}
		return 0;
	}

	bool getValue(vector<int> intRep,T& value)
	{
		int key = LvrHashAlgorithm::DJBHash(intRep);
		typename map<int,T>::iterator it = data.find(key);
		if(it!=data.end())
		{
			value = it->second;
			return true;
		}
		return false;
	}

	void deleteEntry(Atom* atom)
	{
		int key = LvrHashAlgorithm::convertToHash(atom);
		typename map<int,T>::iterator it = data.find(key);
		if(it!=data.end())
		{
			data.erase(it,it++);
		}
	}

	void deleteEntry(vector<int> intRep)
	{
		int key = LvrHashAlgorithm::DJBHash(intRep);
		typename map<int,T>::iterator it = data.find(key);
		if(it!=data.end())
		{
			data.erase(it,it++);
		}
	}

	void clearAll()
	{
		data.clear();
	}

	//KEY RELATED
	/*
	unsigned int DJBHash(vector<int> vals)
	{
		unsigned int hash = 5381;
		for(std::size_t i = 0; i < vals.size(); i++)
		{
			hash = ((hash << 5) + hash) + vals[i];
		}
		return hash;
	}

	int convertToHash(Atom* atom)
	{
		vector<int> signature;
		signature.push_back(atom->symbol->parentId);
		for(unsigned int jj=0;jj<atom->terms.size();jj++)
		{
			if(atom->terms[jj]->domain.size()==1)
			{
				signature.push_back(atom->terms[jj]->domain[0]);
			}
			else
			{
				int sum = atom->terms[jj]->domain.size();
				for(unsigned int kk=0;kk<atom->terms[jj]->domain.size();kk++)
					sum += atom->terms[jj]->domain[kk];
				signature.push_back(sum);
			}
		}
		return DJBHash(signature);
	}
	*/
	int getKey(Atom* atom)
	{
		return LvrHashAlgorithm::convertToHash(atom);
	}

	vector<T> getAllData()
	{
		vector<T> allElems;
		for(typename map<int,T>::iterator it = data.begin();it!=data.end();it++)
		{
			allElems.push_back(it->second);
		}
		return allElems;
	}

	void getAllKeyValuePairs(vector<int>& keys,vector<T>& values)
	{
		for(typename map<int,T>::iterator it = data.begin();it!=data.end();it++)
		{
			keys.push_back(it->first);
			values.push_back(it->second);
		}
	}

	void incrementIntValue(int key)
	{
		typename map<int,T>::iterator it = data.find(key);
		it->second++;
	}

	void setValue(int key,T value)
	{
		typename map<int,T>::iterator it = data.find(key);
		it->second=value;

	}
	
	bool getValue(int key,T& value)
	{
		typename map<int,T>::iterator it = data.find(key);
		if(it == data.end())
			return false;
		value = it->second;
		return true;
	}

	int size()
	{
		return data.size();
	}

	void insert(int key,T value_)
	{
		data.insert(pair<int,T>(key,value_));
	}
};
#endif
