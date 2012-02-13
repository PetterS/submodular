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
			LARGE_INTEGER start_time,stop_time, start_time2, stop_time2;
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

			inline void start2() 
			{ 
				QueryPerformanceCounter(&start_time2);
			}
			double time2()
			{ 
				LARGE_INTEGER freq;
				QueryPerformanceFrequency(&freq);
				return double(stop_time2.QuadPart - start_time2.QuadPart) / double(freq.QuadPart);
			};
			inline double stop2()
			{
				QueryPerformanceCounter(&stop_time2);
				return time2();
			}
		#else
			std::clock_t start_time, stop_time, start_time2, stop_time2;
			inline void start() 
			{ 
				start_time = std::clock();
			}
			double time()
			{ 
				return double(stop_time - start_time) / double(CLOCKS_PER_SEC);
			}
			inline double stop()
			{
				stop_time = std::clock();
				return time();
			}

			inline void start2() 
			{ 
				start_time2 = std::clock();
			}
			double time2()
			{ 
				return double(stop_time2 - start_time2) / double(CLOCKS_PER_SEC);
			}
			inline double stop2()
			{
				stop_time2 = std::clock();
				return time2();
			}
		#endif
	}
}

#endif


