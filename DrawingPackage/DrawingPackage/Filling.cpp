#include "Filling.h"

EdgeRec Filling::getRecord(Point &v1, Point &v2){
	if (v1.y > v2.y)
		std::swap(v1, v2);
	return EdgeRec(v1.x, 1.0*(v1.x - v2.x) / (v1.y - v2.y), v2.y);
}

void Filling::initEdgeRecTable(Point *polygon, int n, std::vector<EdgeRec>  table[]){
	Point v1 = polygon[n - 1];
	for (int i = 0; i < n; ++i){
		Point v2 = polygon[i];
		if (v1.y == v2.y){ v1 = v2; continue; }
		EdgeRec rec = getRecord(v1, v2);
		table[v1.y].push_back(rec);
		v1 = polygon[i];
	}
}

void Filling::polygonFill(HDC hdc, Point *polygon, int n, COLORREF color){
	std::vector<EdgeRec>  table[1000];
	initEdgeRecTable(polygon, n, table);
	int y = 0;
	while (y < 1000 && table[y].size() == 0)++y;
	if (y == 1000)return;
	std::vector<EdgeRec>  ActiveList = table[y];

	while (!ActiveList.empty()){
		std::sort(ActiveList.begin(), ActiveList.end());
		for (int i = 0; i < ActiveList.size() - 1; i += 2){
			int x1 = ceil(ActiveList[i].x);
			int x2 = floor(ActiveList[i + 1].x);
			for (int j = x1; j <= x2; j++){
				SetPixel(hdc, j, y, color);
			}
		}
		y++;
		std::vector<EdgeRec> ::iterator it = ActiveList.begin();
		while (it != ActiveList.end())
		if (y == it->ymax)
			it = ActiveList.erase(it);
		else it++;
		for (std::vector<EdgeRec> ::iterator it = ActiveList.begin(); it != ActiveList.end(); it++)
			it->x += it->minv;
		ActiveList.insert(ActiveList.end(), table[y].begin(), table[y].end());
	}

}

void Filling::initEdgeTable(Edge tbl[]){
	for (int i = 0; i < 1000; ++i){
		tbl[i] = Edge(INT_MAX, INT_MIN);
	}
}

void Filling::constructEdges(Point v1, Point v2, Edge tbl[]){
	if (v1.y > v2.y)
		std::swap(v1, v2);
	if (v1.y == v2.y)return;
	double x = v1.x;
	double minv = 1.0*(v2.x - v1.x) / (v2.y - v1.y);
	for (int y = v1.y; y <= v2.y; ++y){
		tbl[y].xmin = min(tbl[y].xmin, x);
		tbl[y].xmax = max(tbl[y].xmax, x);
		x += minv;
	}
}

void Filling::convexFill(HDC hdc, Point *polygon, int n, COLORREF color){
	Edge tbl[1000];
	initEdgeTable(tbl);
	Point v1 = polygon[n - 1];
	for (int i = 0; i < n; ++i){
		Point v2 = polygon[i];
		constructEdges(v1, v2, tbl);
		v1 = v2;
	}
	for (int y = 0; y < 1000; ++y){
		if (tbl[y].xmin < tbl[y].xmax)
		for (int i = tbl[y].xmin; i <= tbl[y].xmax; ++i){
			SetPixel(hdc, i, y, color);
		}
	}
}


void Filling::FloodFillNonRec(HDC hdc, int x, int y, COLORREF bc, COLORREF fc)
{
	std::stack<Point>st;
	st.push(Point(x, y));
	while (!st.empty())
	{
		Point p = st.top();
		st.pop();
		COLORREF c = GetPixel(hdc, p.x, p.y);
		if (c != bc || c == fc)
		{
			continue;
		}
		SetPixel(hdc, p.x, p.y, fc);
		st.push(Point(p.x + 1, p.y));
		st.push(Point(p.x - 1, p.y));
		st.push(Point(p.x, p.y + 1));
		st.push(Point(p.x, p.y - 1));
	}

}

void Filling::FloodFillRec(HDC hdc, int x, int y, COLORREF bc, COLORREF fc)
{
	COLORREF c = GetPixel(hdc, x, y);
	if (c != bc||c==fc)
	{
		return;
	}
	SetPixel(hdc, x, y, fc);
	FloodFillRec(hdc, x + 1, y, bc, fc);
	FloodFillRec(hdc, x - 1, y, bc, fc);
	FloodFillRec(hdc, x, y + 1, bc, fc);
	FloodFillRec(hdc, x, y - 1, bc, fc);
}