#include "randomgenutil.h"
#include <time.h>
#include <stdlib.h>
LvRandomGenUtil* LvRandomGenUtil::m_pInstance = NULL;

LvRandomGenUtil* LvRandomGenUtil::Instance()
{
	if(!m_pInstance)
	{
		m_pInstance = new LvRandomGenUtil();
	}
	return m_pInstance;
}

LvRandomGenUtil::LvRandomGenUtil()
{
	srand(time(NULL));
	//srand(2350877);
}

int LvRandomGenUtil::getRandomPosition(int maxSize)
{
	return rand()%maxSize;
}

double LvRandomGenUtil::getNormRand()
{
	return (double)rand()/(double)RAND_MAX;
}
