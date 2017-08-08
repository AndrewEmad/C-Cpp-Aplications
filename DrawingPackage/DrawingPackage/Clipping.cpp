#include<Windows.h>
#include<math.h>
#include<algorithm>
#include"Clipping.h"
#include"Line.h"
void Clipping::PointClippingRectangle(HDC hdc,int x, int y, int xmin, int ymin, int xmax, int ymax, COLORREF color)
{
	if (x >= xmin&&x <= xmax&&y >= ymin&&y <= ymax)
	{
		SetPixel(hdc, x, y, color);
	}
}
Point Clipping::VIntersect(double x1, double y1, double x2, double y2, int xEdge)
{
	Point point;
	point.x = xEdge;
	point.y = y1 + (point.x - x1)*(y2 - y1) / (x2 - x1);
	return point;
}
Point Clipping::HTntersect(double x1, double y1, double x2, double y2, int yEdge)
{
	Point point;
	point.y = yEdge;
	point.x = x1 + (point.y - y1)*(x2 - x1) / (y2 - y1);
	return point;
}
OutCode Clipping::getOutCode(double x, double y, int xmin, int ymin, int xmax, int ymax)
{
	OutCode result;
	result.All = 0;
	if (x < xmin)
	{
		result.left = 1;
	}
	else if (x>xmax)
	{
		result.right = 1;
	}
	if (y < ymin)
	{
		result.top = 1;
	}
	else if (y>ymax)
	{
		result.bottom = 1;
	}
	return result;
}
void Clipping::LineClippingRectangle(HDC hdc, double xs, double ys, double xe, double ye, int xmin, int ymin, int xmax, int ymax, COLORREF  color)
{
	OutCode outCode1 = getOutCode(xs, ys, xmin, ymin, xmax, ymax);
	OutCode outCode2 = getOutCode(xe, ye, xmin, ymin, xmax, ymax);
	while ((outCode1.All || outCode2.All) && !(outCode1.All&outCode2.All))
	{
		
		if (outCode1.All)
		{
			Point p;
			if (outCode1.left)
			{
				p = VIntersect(xs, ys, xe, ye, xmin);
			}
			else if (outCode1.right)
			{
				p = VIntersect(xs, ys, xe, ye, xmax);
			}
			else if (outCode1.top)
			{
				p = HTntersect(xs, ys, xe, ye, ymin);
			}
			else if (outCode1.bottom)
			{
				p = HTntersect(xs, ys, xe, ye, ymax);
			}
			xs = p.x;
			ys = p.y; 
			outCode1 = getOutCode(xs, ys, xmin, ymin, xmax, ymax);
		}
		else if (outCode2.All)
		{
			Point p;
			if (outCode2.left)
			{
				p = VIntersect(xs, ys, xe, ye, xmin);
			}
			else if (outCode2.right)
			{
				p = VIntersect(xs, ys, xe, ye, xmax);
			}
			else if (outCode2.top)
			{
				p = HTntersect(xs, ys, xe, ye, ymin);
			}
			else if (outCode2.bottom)
			{
				p = HTntersect(xs, ys, xe, ye, ymax);
			}
			xe = p.x;
			ye = p.y;
			outCode2= getOutCode(xe, ye, xmin, ymin, xmax, ymax);
		}
	}
	if (!outCode1.All&!outCode2.All)
	{
		Line::DrawDDA(hdc, xs, ys, xe, ye, color);
	}
}

void Clipping::PointClippingCircle(HDC hdc, int x, int y, int xc, int yc, double r, COLORREF color)
{
	double result = sqrt((pow((x - xc)*1.0, 2) + pow((y - yc)*1.0, 2)));
	if (result <= r)
	{
		SetPixel(hdc, x, y, color);

	}
}
void Clipping::LineClippingCircle(HDC hdc, double xs, double ys, double xe, double ye, int xc, int yc, double r, COLORREF color)
{
	int dx = xe - xs;
	int dy = ye - ys;
	if (abs(dy) < abs(dx))
	{
		double slope = dy*1.0 / dx;
		if (xs>xe)
		{
			std::swap(xs, xe);
			std::swap(ys, ye);
		}
		int x = xs;
		double y = ys;
		while (x < xe)
		{
			PointClippingCircle(hdc, x, y+0.5, xc, yc, r, color);
			x += 1;
			y += slope;
		}
	}
	else
	{
		double slope = dx*1.0 / dy;
		if (ys>ye)
		{
			std::swap(xs, xe);
			std::swap(ys, ye);
		}
		double x = xs;
		int y = ys;
		while (y < ye)
		{
			PointClippingCircle(hdc,x+0.5,y , xc, yc, r, color);
			x += slope;
			y += 1;
		}
	}
}
