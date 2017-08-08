#include <Windows.h>
#include <math.h>
#include <algorithm>
#include "Circle.h"

void Circle::Draw8Points(HDC hdc, int xc, int yc, int a, int b, COLORREF color) {
	SetPixel(hdc, xc + a, yc + b, color);
	SetPixel(hdc, xc - a, yc + b, color);
	SetPixel(hdc, xc - a, yc - b, color);
	SetPixel(hdc, xc + a, yc - b, color);
	SetPixel(hdc, xc + b, yc + a, color);
	SetPixel(hdc, xc - b, yc + a, color);
	SetPixel(hdc, xc - b, yc - a, color);
	SetPixel(hdc, xc + b, yc - a, color);
}



int Circle::distance(int xs, int ys, int xe, int ye){
	return round(sqrt((xe - xs)*(xe - xs) + (ye - ys)*(ye - ys)));
}


void Circle::DrawParametricCircle(HDC hdc, int xc, int yc, int r, COLORREF c)
{
	double x;
	double y;
	double dceta = 1.0 / r;
	double pi = 22.0 / 7;
	for (double ceta = 0; ceta <= (2 * pi); ceta += dceta)
	{
		x = xc + r*cos(ceta);
		y = yc + r*sin(ceta);
		SetPixel(hdc, round(x), round(y), c);
	}
}


void Circle::DrawDirectCartesian(HDC hdc, int xc, int yc, int r, COLORREF color){
	int y = r;
	for (int x = 0; x <= y; ++x){
		y = round(sqrt(r*r*1.0 - x*x));
		Draw8Points(hdc, xc, yc, x, y, color);
	}

}

void Circle::DrawDirectPolar(HDC hdc, int xc, int yc, int r, COLORREF color){
	int x = 0;
	int y = r;
	double d = 1.0 / r;
	for (double i = 0; !x||x<y; i += d){
		Draw8Points(hdc, xc, yc, x, y, color);
		x = round(r*cos(i));
		y = round(r*sin(i));
	}

}

void Circle::DrawIterativePolar(HDC hdc, int xc, int yc, int r, COLORREF color){
	double x = r;
	double y = 0;
	double dt = 1.0 / r;
	double cdt = cos(dt);
	double sdt = sin(dt);
	while (x >= y){
		Draw8Points(hdc, xc, yc, round(x), round(y), color);
		int tmpX = x;
		x = x*cdt - y*sdt;
		y = tmpX*sdt + y*cdt;
	}
}

void Circle::DrawMidPoint(HDC hdc, int xc, int yc, int r, COLORREF color){
	int y = r;
	int d = 1 - r;
	for (int x = 0; x <= y; ++x){
		Draw8Points(hdc, xc, yc, x, y, color);
		if (d < 0)
			d += 2 * x + 3;
		else{
			y--;
			d += 2 * (x - y) + 5;
		}
	}

}
