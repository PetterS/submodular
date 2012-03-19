/*
*
*   Color
*
*   Allows you to use colors with cout in win32
*   Copyright (c) 2002,2006
*
*   Made by Petter Strandmark 
*
*   Compiles in both Windows and Unix without
*   any changes
*
*/

/*
*  USAGE:
*
*  cout << RED << "This is red." << BLUE << "\nThis is blue.";
*
*
*/


//Retain ANSI/ISO Compability
#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#undef min
#undef max
#endif

#include <iostream>
#include <iomanip>
#include <ctime>

#include "Petter-Color.h"
#include "Petter-Timer.h"


namespace Petter
{
	std::ostream& operator<<(std::ostream& stream,const Color& c)
	{
		stream.flush();
#ifdef _WIN32
		SetConsoleTextAttribute(GetStdHandle(STD_OUTPUT_HANDLE),c.color);
#else
		stream << "\033[0m" << c.str;
#endif
		stream.flush();

		return stream;
	}

#ifdef _WIN32
	const Color NORMAL  = FOREGROUND_RED|FOREGROUND_GREEN|FOREGROUND_BLUE;
	const Color WHITE   = FOREGROUND_RED|FOREGROUND_GREEN|FOREGROUND_BLUE|FOREGROUND_INTENSITY;
	const Color RED     = FOREGROUND_RED|FOREGROUND_INTENSITY;
	const Color DKRED     = FOREGROUND_RED;
	const Color BLUE    = FOREGROUND_BLUE|FOREGROUND_GREEN|FOREGROUND_INTENSITY;
	const Color DKBLUE    = FOREGROUND_BLUE|FOREGROUND_GREEN;
	const Color GREEN   = FOREGROUND_GREEN|FOREGROUND_INTENSITY;
	const Color DKGREEN   = FOREGROUND_GREEN;
	const Color YELLOW  = FOREGROUND_RED|FOREGROUND_GREEN|FOREGROUND_INTENSITY;
	const Color BROWN   = FOREGROUND_RED|FOREGROUND_GREEN;
#else
	const Color NORMAL  = "";
	const Color WHITE   = "\033[37;1m";
	const Color RED     = "\033[31;1m";
	const Color DKRED   = "\033[31m";
	const Color BLUE    = "\033[34;1m";
	const Color DKBLUE  = "\033[34m";
	const Color GREEN   = "\033[32;1m";
	const Color DKGREEN = "\033[32m";
	const Color YELLOW  = "\033[33;1m";
	const Color BROWN   = "\033[33m";
#endif

	std::clock_t start_time = 0;
	bool inside_try_block = false;

	void statusTry(const char* str)
	{
		if (inside_try_block) {
			statusOK();
		}
		std::cerr << std::left << std::setw(40) << str << " [ " << YELLOW << "WAIT" << NORMAL << " ] ";
		inside_try_block = true;
		start2();
	}
	double statusOK()
	{
		stop2();
		if (inside_try_block) {
			double seconds = time2();
			std::cerr << "\b\b\b\b\b\b\b\b" << GREEN << "  OK  " << NORMAL << "]   ";
			std::cerr << /*std::fixed << std::setprecision(1) <<*/ seconds << " s." << std::endl;
			inside_try_block = false;
			return seconds;
		}
		else {
			return 0;
		}
	}
	void statusFailed()
	{
		stop2();
		if (inside_try_block) {
			double seconds = time2();
			std::cerr << "\b\b\b\b\b\b\b\b" << RED << "FAILED" << NORMAL << "]   ";
			std::cerr << /*std::fixed << std::setprecision(1) <<*/ seconds << " s." << std::endl;
			inside_try_block = false;
		}
	}

}

