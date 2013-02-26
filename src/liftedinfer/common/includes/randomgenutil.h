#ifndef __RANDOMGEN_UTIL
#define __RANDOMGEN_UTIL
struct LvRandomGenUtil
{
	static LvRandomGenUtil* Instance();
	int getRandomPosition(int maxSize);
	double getNormRand();
private:
	LvRandomGenUtil();
	static LvRandomGenUtil* m_pInstance;
};
#endif
