/*
 * Chono.h
 *
 *  Created on: 15/set/2016
 *      Author: samuele
 */

#ifndef EVALUATION_CHRONO_H_
#define EVALUATION_CHRONO_H_
#include <unordered_map>
#include <vector>
#include <chrono>
#include <functional>

using namespace std;
using namespace chrono;

class Chrono {
public:
	typedef pair<steady_clock::time_point, steady_clock::time_point> start_end_time;
	Chrono();
	virtual ~Chrono();
	template <typename Obj, typename Function, typename Parameters, typename ReturnType>
	void exe(Obj obj, Function function, Parameters& param, ReturnType& ret, start_end_time& times) {
		times.first = steady_clock::now();
		ret = (obj->*function)(param);
		times.second = steady_clock::now();
	}
	template <typename Obj, typename Function, typename Parameters>
	void exe(Obj obj, Function function, Parameters& param, start_end_time& times) {
		times.first = steady_clock::now();
		(obj->*function)(param);
		times.second = steady_clock::now();
	}
	template <typename Obj, typename Function, typename Parameter1, typename Parameter2>
	void exe1(Obj obj, Function function, Parameter1& param1, Parameter2& param2, start_end_time& times) {
		times.first = steady_clock::now();
		(obj->*function)(param1, param2);
		times.second = steady_clock::now();
	}
};

#endif /* EVALUATION_CHRONO_H_ */
