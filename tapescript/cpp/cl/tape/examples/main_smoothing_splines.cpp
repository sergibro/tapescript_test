//#include <iostream>
//
//#define CL_BASE_SERIALIZER_OPEN
//#define CL_USE_TAPE_SERIALIZATION
//
//#include <cl/tape/impl/detail/enable_ad.hpp>
//#include <cl/tape/tape.hpp>
//#include "impl/utils.hpp"
//#include "impl/polynomial_regression_examples.hpp"
//#include "impl/SmoothingCubicSplines.hpp"
//
//#include "impl/spline_cubic_smooth_slope.hpp"
//#include "impl/smooth_cubic_spline_example.hpp"
//
//#pragma region read_file
//// Reads data from .csv file.
//template<typename Matrix>
//Matrix read(std::string name)
//{
//	typedef typename Matrix::value_type Array;
//	typedef typename Array::value_type T;
//	Matrix values;
//	Array valueline;
//	std::ifstream fin("input/" + name + ".csv");
//	std::string item;
//	for (std::string line; getline(fin, line);)
//	{
//		std::istringstream in(line);
//		while (getline(in, item, ','))
//			valueline.push_back(atof(item.c_str()));
//		values.push_back(valueline);
//		valueline.clear();
//	}
//	return values;
//}
//
//// Reads data from .csv file.
//template <class T>
//inline void read_file(std::string name, std::vector<std::vector<T>>& x, std::vector<T>& y, std::vector<T>& ref)
//{
//	std::ifstream in("input/" + name + ".csv");
//	std::string a;
//	int n = 0;
//	while (!in.eof())
//	{
//		in >> a;
//		n++;
//	}
//	n--;
//	int dim = 0;
//	for (int i = 0; i < a.size(); i++)
//	if (a[i] == 44)
//		dim++;
//	--dim;
//	x.resize(n);
//	for (int i = 0; i < n; i++)
//		x[i].resize(dim);
//	y.resize(n);
//	ref.resize(n);
//	in.close();
//	std::ifstream in1("input/" + name + ".csv");
//	for (int i = 0; i < n; i++)
//	{
//		in1 >> a;
//		std::string value = "";
//		int d = 0;
//		for (int j = 0; j < a.size(); j++)
//		{
//			if (a[j] == 44)
//			{
//				if (d < dim) x[i][d++] = (T)atof(value.c_str());
//				else y[i] = (T)atof(value.c_str());
//				value = "";
//			}
//			else value += a[j];
//			if (j == (a.size() - 1)) ref[i] = (T)atof(value.c_str());
//		}
//	}
//	in1.close();
//}
//#pragma endregion
//
//#pragma region myFuncs
//// Formatted matrix output.
//template<class T>
//void printMatrix(const boost::numeric::ublas::matrix<T> &m, std::ostream & out = std::cout)
//{
//	for (unsigned i = 0; i < m.size1(); ++i)
//	{
//		out << "| ";
//		for (unsigned j = 0; j < m.size2(); ++j)
//		{
//			out << std::setw(14) << std::setprecision(7) << m(i, j) << " | ";
//		}
//		out << "|" << std::endl;
//	}
//	out << std::endl;
//}
//
//template<class T>
//void print(std::vector<T> const& v, std::ostream& out = std::cout)
//{
//	if (v.size() == 0) out << "{}" << std::endl;
//	out.precision(3);
//	out << "{ " << v[0];
//	for (auto i = 1; i < v.size(); i++)
//	{
//		out << ", " << v[i];// << std::endl;
//	}
//	out << " }" << std::endl;
//}
//
//template<class T>
//std::string convertToStr(T *var)
//{
//	std::ostringstream ss;
//	ss << *var;
//	return ss.str();
//}
//
//template<class T>
//void print_csv(std::vector<std::vector<T>> const& xyref, std::vector<T> const& y_res, std::ostream& out = std::cout)
//{
//	for (auto i = 1; i < y_res.size(); ++i)
//		out << xyref.at(i).at(0) << ',' << xyref.at(i).at(1) << ',' << xyref.at(i).at(2) << ',' << y_res.at(i) << std::endl;
//}
//
//template<class T>
//void print_y_res(cl::SmoothingCubicSplines& smoothSpline, std::vector<std::vector<T>> const& for_sort, std::vector<T> const& ref, T const& lambda, std::ostream& out = std::cout)
//{
//	std::vector<T> y_res;
//	smoothSpline.smoothing_spline(lambda);
//	y_res = smoothSpline.get_fitted_value();
//	if (&out != &std::cout) print_csv(for_sort, y_res, out);
//
//	std::cout << lambda << " : " << delta(y_res, ref) << std::endl;//
//}
//
//template<class T>
//T delta(const std::vector<T>& y_res, const std::vector<T>& ref)
//{
//	auto n = ref.size();
//	T d = 0;
//	for (auto i = 0; i < n; ++i)
//	{
//		d += (y_res.at(i) - ref.at(i)) * (y_res.at(i) - ref.at(i));
//	}
//	return sqrt(d / n);
//}
//
//template<class T>
//void test_spline(const bool& write)
//{
//	std::vector<std::string> names1{ "Linear", "LinearClusteredData", "LinearRandomGrid", "PayoffButterfly", "PayoffCallOption", "PayoffCallOptionSell", "ZeroBondAmcExample" },
//		names{ "1_call_geom_brownianGab", "2_RiskReversal", "3_call_geom_brownian", "3_put_geom_brownianGab", "3_Strangle", "9_put_geom_brownian", "11_put_geom_brownian", "17_straddle_geom_brownian", "DownAndInCallBarrierTest1000", "DownAndOutCallBarrierTest1000", "UpAndInCallBarrierTest1000", };
//	for (auto i = 0; i < names.size(); ++i)
//	{
//		std::string name = "/2nd/" + names[i];
//		std::vector<std::vector<T>> xx;
//		std::vector<T> x;
//		std::vector<T> y;
//		std::vector<T> ref;
//		read_file(name, xx, y, ref);
//		for (auto e : xx) x.push_back(e.at(0));
//
//		std::vector<std::vector<T>> for_sort(ref.size());
//		for (auto j = 0; j < ref.size(); ++j)
//		{
//			for_sort.at(j).push_back(x.at(j));
//			for_sort.at(j).push_back(y.at(j));
//			for_sort.at(j).push_back(ref.at(j));
//		}
//		sort(for_sort.begin(), for_sort.end());
//		x.clear();
//		y.clear();
//		ref.clear();
//		for (auto e : for_sort)
//		{
//			x.push_back(e.at(0));
//			y.push_back(e.at(1));
//			ref.push_back(e.at(2));
//		}
//
//		cl::SmoothingCubicSplines smoothSpline;
//		smoothSpline.load(x, y);
//		std::cout << "\tXXX " << names.at(i) << " XXX\nLambda" << std::endl;
//		for (auto j = 0; j < 11; ++j)
//		{
//			double lambda = (j ? (double)j : 0.0001) / 10;
//			std::ofstream out(write ? ("output/2nd(CSV)/" + names.at(i) + "/Lambda_" + convertToStr(&lambda) + "_output.csv") : "output/tmp.txt");
//			print_y_res(smoothSpline, for_sort, ref, lambda, write ? out : std::cout);
//			out.close();
//		}
//	}
//}
//#pragma endregion
//
//int main()
//{
//	auto cl_start = clock();
//#pragma region test_spline
//	//{
//	//	test_spline<double>(false);//true if create outputs, false - just see lambdas
//	//}
//#pragma endregion
//
//#pragma region tmp
//	{
//
//	}
//#pragma endregion
//
//#pragma region testSpline
//	{
//		bool write_base = false;
//		std::ofstream out_base(write_base ? "output/performance_base.csv" : "output/tmp.txt");
//		//cl::spline_cubic_smooth_example(write_base ? out_base : std::cout);
//		out_base.close();
//		
//		bool write_end_value = false;
//		std::ofstream out_end_value(write_end_value ? "output/performance_end_value.csv" : "output/tmp.txt");
//		//cl::spline_cubic_smooth_end_value_example(write_end_value ? out_end_value : std::cout);
//		out_end_value.close();
//
//		bool write_slope = false;
//		std::ofstream out_slope(write_slope ? "output/performance_slope.csv" : "output/tmp.txt");
//		//cl::spline_cubic_smooth_slope_example(write_slope ? out_slope : std::cout);
//		out_slope.close();
//
//		bool write_slope_value = false;
//		std::ofstream out_slope_value(write_slope_value ? "output/performance_slope_value.csv" : "output/tmp.txt");
//		//cl::spline_cubic_smooth_slope_value_example(write_slope_value ? out_slope_value : std::cout);
//		out_slope_value.close();
//
//		bool write_slope_value_diff = true;
//		std::ofstream out_slope_value_diff(write_slope_value_diff ? "output/slope_value_diff.txt" : "output/tmp.txt");
//		cl::spline_cubic_smooth_slope_value_diff_example(write_slope_value_diff ? out_slope_value_diff : std::cout);
//		out_slope_value_diff.close();
//	}
//#pragma endregion
//
//	std::cout << "\nTime: " << clock() - cl_start << " ms" << std::endl;
//	system("pause");
//	return 0;
//}
