#ifndef CURVE_H
#define CURVE_H

#include <Windows.h>
#include <math.h>
#include <algorithm>
#include"Point.h"

#define CURVE_FIRST_DEGREE 14
#define CURVE_SECOND_DEGREE 15
#define CURVE_HERMIT 16
#define CURVE_BEZIER 17
#define CURVE_CARDINAL_SPLINE 18

class Curve{

	public:
		static void DrawFirstDegreeCurve(HDC hdc, Point ps, Point pe, COLORREF c);
		static void DrawSecondDegreeCurve(HDC hdc, Point p1, Point p2, Point p3, COLORREF c);
		static void DrawHermitCurve(HDC hdc, Point p1, Point p2, Point p3, Point p4, COLORREF color);
		static void DrawBezierCurve(HDC hdc, Point p1, Point p2, Point p3, Point p4, COLORREF color);
		static void DrawCardinalSpline(HDC hdc, Point *p, int numOfPoints, COLORREF color);
};

#endif