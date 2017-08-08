#include<vector>
#include<fstream>
#include"Session.h"
#include "Line.h"
#include "Circle.h"
#include "Curve.h"
#include "Filling.h"
#include "Colors.h"
#include "Clipping.h"

using namespace std;

COLORREF Session::BG_Color = C_WHITE;
int Session::BG_Tag = BG_WHITE;
COLORREF Session::FG_Color = C_BLACK;
COLORREF Session::FGC_Color = C_BLACK;
int Session::ClippingAlgorithm = CLIPPING_NONE;
std::vector<Command>Session::Commands;

void Session::setBackColor(COLORREF color){
	BG_Color = color;
}

void Session::setBackTag(int tag){
	BG_Tag = tag;
}

int Session::getBackTag(){
	return BG_Tag;
}

void Session::setClippingAlgorithm(int alg){
	ClippingAlgorithm = alg;
}

int Session::getClippingAlgorithm(){
	return ClippingAlgorithm;
}


void Session::setForeColor(COLORREF color){
	FG_Color = color;
}

COLORREF Session::getBackColor(){
	return BG_Color;
}


void Session::setFGCColor(COLORREF color){
	FGC_Color = color;
}

COLORREF Session::getFGCColor(){
	return FGC_Color;
}


COLORREF Session::getForeColor(){
	return FG_Color;
}

void Session::AddCommand(Command c){
	Commands.push_back(c);
}

void Session::AddCommand(int algorithm ,int numberOfPoints, Point *points, COLORREF color){
	Commands.push_back(Command(algorithm,color, numberOfPoints, points));
}

bool Session::save(string filePath){
	fstream file;
	file.open(filePath, ios::out|ios::binary);
	if (file.is_open()){
		file.write((char*)&BG_Color, sizeof(BG_Color));
		file.write((char*)&BG_Tag, sizeof(BG_Tag));
		file.write((char*)&FGC_Color, sizeof(FGC_Color));
		file.write((char*)&ClippingAlgorithm, sizeof(ClippingAlgorithm));
		int sz = (int)Commands.size();
		file.write((char*)&sz, sizeof(sz));
		for (int i = 0; i < sz; ++i){
			file.write((char*)&Commands[i].Algorithm, sizeof(Commands[i].Algorithm));
			file.write((char*)&Commands[i].Color, sizeof(Commands[i].Color));
			file.write((char*)&Commands[i].numberOfPoints, sizeof(Commands[i].numberOfPoints));
			for (int j = 0; j < Commands[i].numberOfPoints; ++j){
				file.write((char*)&Commands[i].points[j], sizeof(Commands[i].points[j]));
			}
		}
		file.close();
		return true;
	}
	return false;
}

bool Session::load(string filePath){
	fstream file;
	file.open(filePath, ios::out | ios::in | ios::binary);
	if (file.is_open()){
		file.read((char*)&BG_Color, sizeof(BG_Color));
		file.read((char*)&BG_Tag, sizeof(BG_Tag));
		file.read((char*)&FGC_Color, sizeof(FGC_Color));
		file.read((char*)&ClippingAlgorithm, sizeof(ClippingAlgorithm));
		int sz;
		file.read((char*)&sz, sizeof(sz));
		for (int i = 0; i < sz; ++i){
			Command c;
			file.read((char*)&c.Algorithm, sizeof(c.Algorithm));
			file.read((char*)&c.Color, sizeof(c.Color));
			file.read((char*)&c.numberOfPoints, sizeof(c.numberOfPoints));
			c.points = new Point[c.numberOfPoints];
			for (int j = 0; j < c.numberOfPoints; ++j){
				file.read((char*)&c.points[j], sizeof(c.points[j]));
			}
			Commands.push_back(c);
		}
		file.close();
		return true;
	}
	return false;


}

void Session::reDraw(HDC hdc, int xmin, int ymin, int xmax, int ymax, int xc, int yc){

	if (ClippingAlgorithm == CLIPPING_NONE){
		Line::DrawDirect(hdc, xmin, ymin, xmin, ymax, Session::getBackColor());
		Line::DrawDirect(hdc, xmin, ymin, xmax, ymin, Session::getBackColor());
		Line::DrawDirect(hdc, xmax, ymin, xmax, ymax, Session::getBackColor());
		Line::DrawDirect(hdc, xmin, ymax, xmax, ymax, Session::getBackColor());
		Circle::DrawMidPoint(hdc, xc, yc, (xmax - xmin) / 2, Session::getBackColor());
	}
	else if (ClippingAlgorithm == CLIPPING_CIRCLE){
		Line::DrawDirect(hdc, xmin, ymin, xmin, ymax, Session::getBackColor());
		Line::DrawDirect(hdc, xmin, ymin, xmax, ymin, Session::getBackColor());
		Line::DrawDirect(hdc, xmax, ymin, xmax, ymax, Session::getBackColor());
		Line::DrawDirect(hdc, xmin, ymax, xmax, ymax, Session::getBackColor());
		Circle::DrawMidPoint(hdc, xc, yc, (xmax - xmin) / 2, Session::getFGCColor());
	}
	else if (ClippingAlgorithm == CLIPPING_RECTANGLE){
		Circle::DrawMidPoint(hdc, xc, yc, (xmax - xmin) / 2, Session::getBackColor());
		Line::DrawDirect(hdc, xmin, ymin, xmin, ymax, Session::getFGCColor());
		Line::DrawDirect(hdc, xmin, ymin, xmax, ymin, Session::getFGCColor());
		Line::DrawDirect(hdc, xmax, ymin, xmax, ymax, Session::getFGCColor());
		Line::DrawDirect(hdc, xmin, ymax, xmax, ymax, Session::getFGCColor());
	}
	
	for (int i = 0; i < Commands.size(); ++i){

		Command command = Commands[i];

		switch (command.Algorithm){
			case D_POINT:
				if (ClippingAlgorithm == CLIPPING_NONE){
					SetPixel(hdc, command.points[0].x, command.points[0].y, command.Color);
				}
				else if (ClippingAlgorithm = CLIPPING_RECTANGLE){
					Clipping::PointClippingRectangle(hdc, command.points[0].x, command.points[0].y, xmin, ymin, xmax, ymax, command.Color);
				}
				else if (ClippingAlgorithm = CLIPPING_CIRCLE){
					Clipping::PointClippingCircle(hdc, command.points[0].x, command.points[0].y, xc, yc, (xmax - xmin) / 2, command.Color);
				}
				break;


			case LINE_PARAMETRIC:
				if (ClippingAlgorithm == CLIPPING_NONE)
					Line::DrawParametric(hdc, command.points[0].x, command.points[0].y, command.points[1].x, command.points[1].y, command.Color);
				else if (ClippingAlgorithm == CLIPPING_RECTANGLE)
					Clipping::LineClippingRectangle(hdc, command.points[0].x, command.points[0].y, command.points[1].x, command.points[1].y, xmin, ymin, xmax, ymax, command.Color);
				else if (ClippingAlgorithm == CLIPPING_CIRCLE)
					Clipping::LineClippingCircle(hdc, command.points[0].x, command.points[0].y, command.points[1].x, command.points[1].y, xc, yc, (xmax - xmin) / 2, command.Color);

					break;

			case LINE_DIRECT:
				Line::DrawDirect(hdc, command.points[0].x, command.points[0].y, command.points[1].x, command.points[1].y, command.Color);
				break;

			case LINE_DDA:
				Line::DrawDDA(hdc, command.points[0].x, command.points[0].y, command.points[1].x, command.points[1].y, command.Color);

				break;

			case LINE_MIDPOINT:
				Line::DrawMidPoint(hdc, command.points[0].x, command.points[0].y, command.points[1].x, command.points[1].y, command.Color);
				break;

			case CIRCLE_PARAMETRIC:
				Circle::DrawParametricCircle(hdc, command.points[0].x, command.points[0].y, Circle::distance(command.points[0].x, command.points[0].y, command.points[1].x, command.points[1].y), command.Color);

				break;

			case CIRCLE_DIRECT_CARTESIAN:
				Circle::DrawDirectCartesian(hdc, command.points[0].x, command.points[0].y, Circle::distance(command.points[0].x, command.points[0].y, command.points[1].x, command.points[1].y), command.Color);

				break;

			case CIRCLE_DIRECT_POLAR:
				Circle::DrawDirectPolar(hdc, command.points[0].x, command.points[0].y, Circle::distance(command.points[0].x, command.points[0].y, command.points[1].x, command.points[1].y), command.Color);
				break;

			case CIRCLE_ITERATIVE_POLAR:
				Circle::DrawIterativePolar(hdc, command.points[0].x, command.points[0].y, Circle::distance(command.points[0].x, command.points[0].y, command.points[1].x, command.points[1].y), command.Color);
				break;

			case CIRCLE_MIDPOINT:
				Circle::DrawMidPoint(hdc, command.points[0].x, command.points[0].y, Circle::distance(command.points[0].x, command.points[0].y, command.points[1].x, command.points[1].y), command.Color);

				break;

			case CURVE_FIRST_DEGREE:
				Curve::DrawFirstDegreeCurve(hdc, command.points[0], command.points[1], command.Color);

				break;

			case CURVE_SECOND_DEGREE:
				Curve::DrawSecondDegreeCurve(hdc, command.points[0], command.points[1], command.points[2], command.Color);

				break;

			case CURVE_HERMIT:
				Curve::DrawHermitCurve(hdc, command.points[0], command.points[1], command.points[2], command.points[3], command.Color);

				break;

			case CURVE_BEZIER:
				Curve::DrawBezierCurve(hdc, command.points[0], command.points[1], command.points[2], command.points[3], command.Color);

				break;

			case CURVE_CARDINAL_SPLINE:
				Curve::DrawCardinalSpline(hdc, command.points, command.numberOfPoints, command.Color);

				break;

			case FILLED_SHAPES_CONVEX:
				Filling::convexFill(hdc, command.points, command.numberOfPoints, command.Color);

				break;

			case FILLED_SHAPES_POLYGON:
				Filling::polygonFill(hdc, command.points, command.numberOfPoints, command.Color);
				break;


			case FILLING_FLOODFILL_REC:
				Filling::FloodFillRec(hdc, command.points[0].x, command.points[0].y, Session::getBackColor(), Session::getForeColor());
				break;

			case FILLING_FLOODFILL_NON_REC:
				Filling::FloodFillNonRec(hdc, command.points[0].x, command.points[0].y, Session::getBackColor(), Session::getForeColor());
				break;
		}
	}
}


void Session::clear(){
	Commands.clear();
}