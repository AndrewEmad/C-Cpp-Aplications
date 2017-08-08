#include <Windows.h>
#include <math.h>
#include <algorithm>
#include "Curve.h"


void Curve::DrawFirstDegreeCurve(HDC hdc, Point ps, Point pe, COLORREF c)
{
	int xCofficient[2] = { pe.x - ps.x, ps.x };
	int yCofficient[2] = { pe.y - ps.y, ps.y };
	double x, y;
	double dt = (max(abs(pe.x - ps.x), abs(pe.y - ps.y)));
	dt = 1.0 / dt;
	for (double t = 0; t < 1; t += dt)
	{
		x = xCofficient[0] * t + xCofficient[1];
		y = yCofficient[0] * t + yCofficient[1];
		SetPixel(hdc, round(x), round(y), c);
	}
}

void Curve::DrawSecondDegreeCurve(HDC hdc, Point p1, Point p2, Point p3, COLORREF c)
{
	int xCofficient[3] = { 2 * p3.x - 4 * p2.x + 2 * p1.x, 4 * p2.x - p3.x - 3 * p1.x, p1.x };
	int yCofficient[3] = { 2 * p3.y - 4 * p2.y + 2 * p1.y, 4 * p2.y - p3.y - 3 * p1.y, p1.y };
	double dt = (max(abs(p1.x - p2.x), abs(p1.y - p2.y)));
	dt = max(dt, (abs(p2.x - p3.x), abs(p2.y - p3.y)));
	dt = max(dt, (abs(p1.x - p3.x), abs(p1.y - p3.y)));
	dt = 1.0 / dt;
	double x, y;
	for (double t = 0; t < 1; t += 0.0001)
	{
		x = xCofficient[0] * pow(t, 2) + xCofficient[1] * t + xCofficient[2];
		y = yCofficient[0] * pow(t, 2) + yCofficient[1] * t + yCofficient[2];
		SetPixel(hdc, round(x), round(y), c);
	}

}

void Curve::DrawHermitCurve(HDC hdc, Point p1, Point p2, Point p3, Point p4, COLORREF color)
{

	int xcoff[4] = { 2 * p1.x + p2.x - 2 * p3.x + p4.x, -3 * p1.x - 2 * p2.x + 3 * p3.x - p4.x, p2.x, p1.x };
	int ycoff[4] = { 2 * p1.y + p2.y - 2 * p3.y + p4.y, -3 * p1.y - 2 * p2.y + 3 * p3.y - p4.y, p2.y, p1.y };

	for (double t = 0; t < 1; t += 0.001){
		int x = xcoff[0] * t*t*t + xcoff[1] * t*t + xcoff[2] * t + xcoff[3] + 0.5;
		int y = ycoff[0] * t*t*t + ycoff[1] * t*t + ycoff[2] * t + ycoff[3] + 0.5;
		SetPixel(hdc, x, y, color);
	}
}

void Curve::DrawBezierCurve(HDC hdc, Point p1, Point p2, Point p3, Point p4, COLORREF color)
{
	Point T1, T2;
	T1.x = 3 * (p2.x - p1.x);
	T1.y = 3 * (p2.y - p1.y);
	T2.x = 3 * (p4.x - p3.x);
	T2.y = 3 * (p4.y - p3.y);
	DrawHermitCurve(hdc, p1, T1, p3, T2, color);
}


void Curve::DrawCardinalSpline(HDC hdc, Point *p, int numOfPoints, COLORREF color)
{
	double c1 = 1 - 0.001;
	Point T0, T1;
	T0.x = c1*(p[2].x - p[0].x);
	T0.y = c1*(p[2].y - p[0].y);
	for (int i = 2; i < numOfPoints - 1; ++i){
		T1.x = c1*(p[i+1].x - p[i-1].x);
		T1.y = c1*(p[i+1].y - p[i-1].y);
		DrawHermitCurve(hdc,p[i-1],T0,p[i],T1,color);
		T0 = T1;
	}
}