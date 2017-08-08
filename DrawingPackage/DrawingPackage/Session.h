#ifndef SESSION_H
#define SESSION_H
#include<vector>
#include<Windows.h>
#include"Point.h"

struct Command{
	int Algorithm;
	COLORREF Color;
	int numberOfPoints;
	Point *points;
	Command(){}
	Command(int algorithm, COLORREF clr, int nop, Point*p){
		Algorithm = algorithm;
		Color = clr;
		numberOfPoints = nop;
		points = new Point[numberOfPoints];
		for (int i = 0; i < numberOfPoints; ++i){
			points[i].x = p[i].x;
			points[i].y = p[i].y;
		}
	}
};

class Session{
	private:
		static COLORREF BG_Color;
		static int BG_Tag;
		static COLORREF FG_Color;
		static COLORREF FGC_Color;
		static int ClippingAlgorithm;
		static std::vector<Command>Commands;
	public:
		static void setBackColor(COLORREF color);
		static void setBackTag(int tag);
		static void setClippingAlgorithm(int alg);
		static void setForeColor(COLORREF color);
		static void setFGCColor(COLORREF color);
		static COLORREF getFGCColor();
		static COLORREF getBackColor();
		static int getBackTag();
		static int getClippingAlgorithm();
		static COLORREF getForeColor();
		static void AddCommand(Command c);
		static void AddCommand(int algorithm, int numberOfPoints, Point *points, COLORREF color);
		static bool save(std::string filePath);
		static bool load(std::string filePath);
		static void reDraw(HDC hdc, int xmin, int ymin, int xmax, int ymax, int xc, int yc);
		static void clear();
};






#endif