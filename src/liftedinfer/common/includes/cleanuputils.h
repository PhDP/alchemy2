#ifndef __CLEANUP_UTILS_H
#define __CLEANUP_UTILS_H

template<typename Container>
void cleanup(Container& c)
{
	while(!c.empty())
	{
		
		if(c.back())
		{
			delete c.back();
			c.back() = NULL;
		}
		c.pop_back();
	}
}

template<typename Container>
void removeItem(Container& c, int itemIndex)
{
	
	if(c.at(itemIndex))
	{
		delete c.at(itemIndex);
		c.at(itemIndex) = NULL;
	}
	c.erase(c.begin()+itemIndex);
}


template<typename Container>
void removeUnusedItems(Container& oldC, Container& newC)
{
	while(!oldC.empty())
	{
		if(!oldC.back())
			continue;
		bool found = false;
		for(unsigned int i=0;i<newC.size();i++)
		{
			if(oldC.back() == newC[i])
			{
				found =true;
				break;
			}
		}
		if(!found)
		{
			//element does not occur in new container, delete it
			delete oldC.back();
			oldC.back() = NULL;

		}
	oldC.pop_back();
	}
}

template<typename Container>
void cleanupItems(Container& c,int startIndex = -1,int endIndex = -1)
{
	int start;
	startIndex == -1?start = 0:start = startIndex;
	int end;
	endIndex == -1?end = c.size():end = endIndex;
	for(unsigned int m=start;m<end;m++)
		cleanup(c[m]);
}

template<typename Container>
void removeItems(Container& c,int startIndex = -1,int endIndex = -1)
{
	int start;
	startIndex == -1?start = 0:start = startIndex;
	int end;
	endIndex == -1?end = c.size():end = endIndex;
	for(unsigned int m=start;m<end;m++)
	{
		delete c[m];
	}
	c.erase(c.begin()+start,c.begin()+end);
}
#endif