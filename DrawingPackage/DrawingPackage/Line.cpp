#include <Windows.h>
#include <math.h>
#include <algorithm>
#include "Line.h"

void Line::DrawParametric(HDC hdc, int xs, int ys, int xe, int ye, COLORREF color){
	int n = max(abs(xe - xs), abs(ye - ys));
	double dt = 1.0 / n;
	if (n == abs(xe - xs)){
		if (xe < xs)
		{
			std::swap(xe, xs);
			std::swap(ye, ys);
		}
		double dy = (ye - ys)*1.0 / n;
		double y = ys;
		for (int x = xs; x <= xe; ++x, y += dy){
			SetPixel(hdc, x, round(y), color);
		}


	}
	else{
		if (ye < ys)
		{
			std::swap(xe, xs);
			std::swap(ye, ys);
		}
		double dx = (xe - xs)*1.0 / n;
		double x = xs;
		for (int y = ys; y <= ye; ++y, x += dx){
			SetPixel(hdc, round(x), y, color);
		}


	}
}

void Line::DrawDirect(HDC hdc, int xs, int ys, int xe, int ye, COLORREF color){
	int dx = xe - xs;
	int dy = ye - ys;
	if (abs(dx) > abs(dy)){
		if (xe < xs){
			std::swap(xe, xs);
			std::swap(ye, ys);
		}
		double m = dy*1.0 / dx;

		for (int x = xs; x < xe; ++x){
			SetPixel(hdc, x, round(ys + m*(x - xs)), color);
		}
	}

	else{
		if (ye < ys){
			std::swap(xe, xs);
			std::swap(ye, ys);
		}
		double m = dx*1.0 / dy;

		for (int y = ys; y < ye; ++y){
			SetPixel(hdc, round(xs + m*(y - ys)), y, color);
		}
	}
}

void Line::DrawDDA(HDC hdc, int xs, int ys, int xe, int ye, COLORREF color){
	int dx = xe - xs;
	int dy = ye - ys;
	if (abs(dx) > abs(dy)){
		if (xe < xs){
			std::swap(xe, xs);
			std::swap(ye, ys);
		}
		double m = dy*1.0 / dx;
		double y = ys;
		for (int x = xs; x < xe; ++x, y += m){
			SetPixel(hdc, x, round(y), color);
		}
	}

	else{
		if (ye < ys){
			std::swap(xe, xs);
			std::swap(ye, ys);
		}
		double m = dx*1.0 / dy;
		double x = xs;
		for (int y = ys; y < ye; ++y, x += m){
			SetPixel(hdc, round(x), y, color);

		}
	}
}

void Line::DrawMidPoint(HDC hdc, int xs, int ys, int xe, int ye, COLORREF color){
	int dx = xe - xs;
	int dy = ye - ys;
	if (abs(dx) > abs(dy)){
		if (xe < xs){
			std::swap(xe, xs);
			std::swap(ye, ys);
			dx = xe - xs;
			dy = ye - ys;
		}
		if (ys < ye){
			int d = dx - 2 * dy;
			int d1 = 2 * (dx - dy);
			int d2 = -2 * dy;

			for (int x = xs; x <= xe; ++x){
				SetPixel(hdc, x, ys, color);
				if (d < 0){
					d += d1;
					ys++;
				}
				else
					d += d2;
			}
		}
		else{
			int d = -dx - 2 * dy;
			int d1 = -2 * (dx + dy);
			int d2 = -2 * dy;

			for (int x = xs; x <= xe; ++x){
				SetPixel(hdc, x, ys, color);
				if (d > 0){
					d += d1;
					ys--;
				}
				else
					d += d2;
			}

		}
	}

	else{
		if (ye < ys){
			std::swap(xe, xs);
			std::swap(ye, ys);
			dx = xe - xs;
			dy = ye - ys;
		}

		if (xs < xe){

			int d = 2 * dx - dy;
			int d1 = 2 * (dx - dy);
			int d2 = 2 * dx;

			for (int y = ys; y <= ye; ++y){
				SetPixel(hdc, xs, y, color);
				if (d > 0){
					d += d1;
					xs++;
				}
				else
					d += d2;
			}
		}

		else{
			int d = 2 * dx + dy;
			int d1 = 2 * (dx + dy);
			int d2 = 2 * dx;

			for (int y = ys; y <= ye; ++y){
				SetPixel(hdc, xs, y, color);
				if (d < 0){
					d += d1;
					xs--;
				}
				else
					d += d2;
			}
		}
	}
}

