#ifndef CIRCLE_H
#define CIRCLE_H

#include <Windows.h>
#include <math.h>
#include <algorithm>

#define CIRCLE_PARAMETRIC 9
#define CIRCLE_DIRECT_CARTESIAN 10
#define CIRCLE_DIRECT_POLAR 11
#define CIRCLE_ITERATIVE_POLAR 12
#define CIRCLE_MIDPOINT 13

class Circle{
	private:
		static void Draw8Points(HDC hdc, int xc, int yc, int a, int b, COLORREF color);



	public:
		static void DrawParametricCircle(HDC hdc, int xc, int yc, int r, COLORREF c);

		static void DrawDirectCartesian(HDC hdc, int xc, int yc, int r, COLORREF color);

		static void DrawDirectPolar(HDC hdc, int xc, int yc, int r, COLORREF color);

		static void DrawIterativePolar(HDC hdc, int xc, int yc, int r, COLORREF color);

		static void DrawMidPoint(HDC hdc, int xc, int yc, int r, COLORREF color);

		static int distance(int xs, int ys, int xe, int ye);


};

#endif