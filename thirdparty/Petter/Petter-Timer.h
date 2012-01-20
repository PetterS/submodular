#ifndef PETTER_TIMER_H
#define PETTER_TIMER_H

#ifdef _WIN32
	#define WIN32_LEAN_AND_MEAN
	#include <windows.h>
#else
	#include <ctime>
#endif

namespace Petter 
{
	namespace
	{
		#ifdef _WIN32
			LARGE_INTEGER start_time,stop_time;
			inline void start() 
			{ 
				QueryPerformanceCounter(&start_time);
			}
			double time()
			{ 
				LARGE_INTEGER freq;
				QueryPerformanceFrequency(&freq);
				return double(stop_time.QuadPart - start_time.QuadPart) / double(freq.QuadPart);
			};
			inline double stop()
			{
				QueryPerformanceCounter(&stop_time);
				return time();
			}
		#else
			std::clock_t start_time, stop_time;
			inline void start() 
			{ 
				start_time = std::clock();
			}
			double time()
			{ 
				return double(stop_time - start_time) / double(CLOCKS_PER_SEC);
			};
			inline double stop()
			{
				stop_time = std::clock();
				return time();
			}
		#endif
	}
}

#endif


