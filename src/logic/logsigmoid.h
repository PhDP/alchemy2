#ifndef LOG_SIGMOID_DEC_2007
#define LOG_SIGMOID_DEC_2007


//currently LogSigMoid(t, a) is taking the form -log(1+exp(a*(t-x)))
class LogSigMoid
{
	public:
		LogSigMoid()
		{
			strictness_ = 0;
			threshold_ = 0;
		}

		LogSigMoid(double strictness, double threshold): strictness_(strictness), threshold_(threshold)
		{}

		double funcValue(double x) //compute function value
		{
			return -log(1+exp(strictness_*(threshold_-x)));
		}

		double firstDerivative(double x) //compute first derivative
		{
			// first derivative is of the form 
			// a*exp(a*(t-x))/(1+exp(a*(t-x)))
			return strictness_*exp(strictness_*(threshold_-x))/(1+exp(strictness_*(threshold_-x)));
		}

		double solveConstraint(double th) // solve constraint of the form logSigmoid(a, t, x) > th
		{
			if (th >=0 ) return NOSOL;

			double root = threshold_ - log(exp(-1*th) - 1)/strictness_;

			return root;
		}
	public:
		double strictness_; // x > t if strictness_ > 0, x < t if strictness < 0
		double threshold_; //
};
#endif
