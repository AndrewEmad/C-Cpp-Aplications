#ifndef LINE_H
#define LINE_H

#include <Windows.h>
#include <math.h>
#include <algorithm>


#define LINE_PARAMETRIC 5
#define LINE_DIRECT 6
#define LINE_DDA 7
#define LINE_MIDPOINT 8

class Line{
	
	public:
		static void DrawParametric(HDC hdc, int xs, int ys, int xe, int ye, COLORREF color);
		static void DrawDirect(HDC hdc, int xs, int ys, int xe, int ye, COLORREF color);
		static void DrawDDA(HDC hdc, int xs, int ys, int xe, int ye, COLORREF color);
		static void DrawMidPoint(HDC hdc, int xs, int ys, int xe, int ye, COLORREF color);

};

#endif