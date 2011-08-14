/*
*
*   Color
*
*   Allows you to use colors with cout in win32
*   Copyright (c) 2002,2006,2010
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

#ifndef COLOR_PETTER_H
#define COLOR_PETTER_H


//Retain ANSI/ISO Compability
#ifdef _WIN32
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#undef min
#undef max
#endif

#include <iostream>
#include <ctime>


namespace Petter
{
	class Color
	{
		friend std::ostream& operator<<(std::ostream& stream,const Color& c);

	public:
#ifdef _WIN32
		Color(unsigned short c): color(c) {}      
		unsigned short color;    
#else
		Color(const char* s): str(s) {}
		const char* str;
#endif   
	};

	std::ostream& operator<<(std::ostream& stream,const Color& c);

	extern const Color NORMAL;
	extern const Color WHITE;
	extern const Color RED;
	extern const Color DKRED;
	extern const Color BLUE;
	extern const Color DKBLUE;
	extern const Color GREEN;
	extern const Color DKGREEN;
	extern const Color YELLOW;
	extern const Color BROWN;

	extern std::clock_t start_time;
	extern bool inside_try_block;

	void statusTry(const char* str);
	double statusOK();
	void statusFailed();
}



#endif //ifndef 

