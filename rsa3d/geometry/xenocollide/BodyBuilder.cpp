/*
XenoCollide Collision Detection and Physics Library
Copyright (c) 2007-2014 Gary Snethen http://xenocollide.com

This software is provided 'as-is', without any express or implied warranty.
In no event will the authors be held liable for any damages arising
from the use of this software.
Permission is granted to anyone to use this software for any purpose,
including commercial applications, and to alter it and redistribute it freely,
subject to the following restrictions:

1. The origin of this software must not be misrepresented; you must
not claim that you wrote the original software. If you use this
software in a product, an acknowledgment in the product documentation
would be appreciated but is not required.
2. Altered source versions must be plainly marked as such, and must
not be misrepresented as being the original software.
3. This notice may not be removed or altered from any source distribution.
*/



#include "BodyBuilder.h"
#include "MapPtr.h"
#include <list>

//////////////////////////////////////////////////////////////

BodyBuilder::BodyBuilder()
{

}

//////////////////////////////////////////////////////////////

double BodyBuilder::getMaxRadius()
{
	XCShape& s = *mShapeStack.back();
	CollideGeometry *collideModel = s.geom;
	double radiusNegX = std::abs(collideModel->GetSupportPoint( Vector<3>({-1, 0, 0}) )[0]);
	double radiusPosX = std::abs(collideModel->GetSupportPoint( Vector<3>({1, 0, 0}) )[0]);
	double radiusNegY = std::abs(collideModel->GetSupportPoint( Vector<3>({0, -1, 0}) )[1]);
	double radiusPosY = std::abs(collideModel->GetSupportPoint( Vector<3>({0,  1, 0}) )[1]);
	double radiusNegZ = std::abs(collideModel->GetSupportPoint( Vector<3>({0, 0, -1}) )[2]);
	double radiusPosZ = std::abs(collideModel->GetSupportPoint( Vector<3>({0, 0,  1}) )[2]);

	Vector<3> maxRadiusVector
	({
		std::max(radiusNegX, radiusPosX),
		std::max(radiusNegY, radiusPosY),
		std::max(radiusNegZ, radiusPosZ)
	});

	return maxRadiusVector.norm();
}

MapPtr<CollideGeometry> BodyBuilder::getCollideGeometry(){
    ValidateMsg(mShapeStack.size()>=1, "Empty shape stack");
	XCShape& s = *mShapeStack.back();
	return s.geom;
}


BodyBuilder::~BodyBuilder()
{
	mShapeStack.clear();
}

void BodyBuilder::axis(double x){
	MapPtr<CollideGeometry> geom = new CollideSegment(x);
	mShapeStack.push_back( new XCShape(geom) );
	mShapeStack.push_back( new XCShape(geom, Quat(0, 0, std::sin(M_PI/4), std::cos(M_PI/4)), Vector<3>({0, 0, 0}) ) );
	mShapeStack.push_back( new XCShape(geom, Quat(0, std::sin(M_PI/4), 0, std::cos(M_PI/4)), Vector<3>({0, 0, 0}) ) );
}

void BodyBuilder::box(double x, double y, double z){
	MapPtr<CollideGeometry> geom = new CollideBox(Vector<3>({x, y, z}));
	mShapeStack.push_back( new XCShape(geom) );
}

void BodyBuilder::clear(){
	this->mShapeStack.clear();
}

void BodyBuilder::diff(){
	ValidateMsg(mShapeStack.size() >= 2, "Insufficient shapes number on stack. At least two required");
	MapPtr<XCShape> shape2 = mShapeStack.back();
	mShapeStack.pop_back();

	MapPtr<XCShape> shape1 = mShapeStack.back();
	mShapeStack.pop_back();

	MapPtr<XCShape> newShape = new XCShape();
	newShape->geom = new CollideDiff(shape1->geom, shape1->q, shape1->x, shape2->geom, shape2->q, shape2->x);
	mShapeStack.push_back(newShape);
}

void BodyBuilder::disc(double x){
	MapPtr<CollideGeometry> geom = new CollideDisc(x);
	mShapeStack.push_back( new XCShape(geom) );
}

void BodyBuilder::disc(double x, double y){
	MapPtr<CollideGeometry> geom = new CollideEllipse(x, y);
	mShapeStack.push_back( new XCShape(geom) );
}

void BodyBuilder::dup(size_t n){
	ValidateMsg(mShapeStack.size() >= n, "Insufficient shapes number on stack. At least n required");
	std::list< MapPtr<XCShape> >::iterator it = mShapeStack.end();
	for (size_t i=0; i < n; i++){
			it--;
	}
	for (size_t i=0; i < n; i++){
		MapPtr<XCShape> s = new XCShape();
		*s = **it;
		mShapeStack.push_back( s );
		it++;
	}
}

void BodyBuilder::football(double l, double w){
	MapPtr<CollideGeometry> geom = new CollideFootball(l,w);
	mShapeStack.push_back( new XCShape(geom) );
}

void BodyBuilder::move(double x, double y, double z){
	ValidateMsg(mShapeStack.size() > 0, "Empty stack. No shape to move");
	mShapeStack.back()->x += Vector<3>({x, y, z});
}

void BodyBuilder::point(double x, double y, double z){
	MapPtr<CollideGeometry> geom = new CollidePoint(Vector<3>({x, y, z}));
	mShapeStack.push_back( new XCShape(geom) );
}

void BodyBuilder::pop(){
	ValidateMsg(mShapeStack.size() > 0, "Empty stack. No shape to pop");
	mShapeStack.pop_back();
}

void::BodyBuilder::rect(double x, double y){
	MapPtr<CollideGeometry> geom = new CollideRectangle(x, y);
	mShapeStack.push_back( new XCShape(geom) );
}

void BodyBuilder::rot(double x, double y, double z){
	ValidateMsg(mShapeStack.size() > 0, "Empty stack. No shape to rotate");
	Quat quatLog( x*M_PI/180.0, y*M_PI/180.0, z*M_PI/180.0, 0);
	Quat q = Quat::QuatExp(quatLog);
	mShapeStack.back()->q = q * mShapeStack.back()->q;
}

void BodyBuilder::saucer(double r, double t){
	MapPtr<CollideGeometry> geom = new CollideSaucer(r,t);
	mShapeStack.push_back( new XCShape(geom) );
}

void BodyBuilder::segment(double l){
	MapPtr<CollideGeometry> geom = new CollideSegment(l);
	mShapeStack.push_back( new XCShape(geom) );
}

void BodyBuilder::sphere(double r){
	MapPtr<CollideGeometry> geom = new CollideSphere(r);
	mShapeStack.push_back( new XCShape(geom) );
}

void BodyBuilder::sphere(double rx, double ry, double rz){
	MapPtr<CollideGeometry> geom = new CollideEllipsoid(Vector<3>({rx, ry, rz}));
	mShapeStack.push_back( new XCShape(geom) );
}

void BodyBuilder::swap(){
	ValidateMsg(mShapeStack.size() >= 2, "Insufficient shapes number on stack. At least two required");
	MapPtr<XCShape> s1 = mShapeStack.back();
	mShapeStack.pop_back();

	MapPtr<XCShape> s2 = mShapeStack.back();
	mShapeStack.pop_back();

	mShapeStack.push_back(s1);
	mShapeStack.push_back(s2);
}

void BodyBuilder::sweep(){
	ValidateMsg(mShapeStack.size() >= 2, "Insufficient shapes number on stack. At least two required");
	MapPtr<XCShape> shape1 = mShapeStack.back();
	mShapeStack.pop_back();

	MapPtr<XCShape> shape2 = mShapeStack.back();
	mShapeStack.pop_back();

	MapPtr<XCShape> newShape = new XCShape();
	newShape->geom = new CollideSum(shape1->geom, shape1->q, shape1->x, shape2->geom, shape2->q, shape2->x);
	mShapeStack.push_back(newShape);
}

void BodyBuilder::wrap(){
	ValidateMsg(mShapeStack.size() >= 2, "Insufficient shapes number on stack. At least two required");
	MapPtr<XCShape> shape1 = mShapeStack.back();
	mShapeStack.pop_back();

	MapPtr<XCShape> shape2 = mShapeStack.back();
	mShapeStack.pop_back();

	MapPtr<XCShape> newShape = new XCShape();
	newShape->geom = new CollideMax(shape1->geom, shape1->q, shape1->x, shape2->geom, shape2->q, shape2->x);
	mShapeStack.push_back(newShape);
}

//////////////////////////////////////////////////////////////


void BodyBuilder::ProcessCommand(std::string& cmd)
{
	for (size_t i=0; i < cmd.size(); i++){
		cmd[i] = std::tolower(cmd[i]);
	}

	if (cmd == "axis"){
		double x = 10;
		sscanf(cmd.c_str(), " axis %lf", &x);
		this->axis(x);
		Status("Axis created.");
	}
	else if (cmd.substr(0,3) == "box"){
		double x = 5;
		double y = 5;
		double z = 5;
		sscanf(cmd.c_str(), "box %lf %lf %lf", &x, &y, &z);
		this->box(x, y, z);
		Status("Box created.");
	}
	else if (cmd == "diff"){
		if (mShapeStack.size() >= 2){
			this->diff();
			Status("Minkowski difference created.");
		}else{
			Status("Minkowski difference requires two shapes.");
		}
	}
	else if (cmd.substr(0,4) == "disc"){
		double x = 5;
		double y = 5;
		double count = sscanf(cmd.c_str(), "disc %lf %lf", &x, &y);
		if (count == 2){
			this->disc(x, y);
		}else{
			this->disc(x);
		}
		Status("Disc created.");
	}
	else if (cmd.substr(0, 3) == "dup"){
		size_t n = 1;
		sscanf(cmd.c_str(), "dup %lu", &n);
		if (mShapeStack.size() >= n){
			this->dup(n);
			Status("Shape(s) duplicated.");
		}else{
			Status("Too few shapes on stack.");
		}
	}
	else if (cmd.substr(0, 7) == "execute" || cmd.substr(0, 3) == "run")
	{
		char name[100];
		if (sscanf(cmd.c_str(), " execute %s", name) == 1 || sscanf(cmd.c_str(), " run %s", name) == 1)
		{
			std::string filename = std::string("scripts\\") + name;
			FILE* fp = fopen(filename.c_str(), "rt");
			if (!fp)
			{
				filename += ".txt";
				fp = fopen(filename.c_str(), "rt");
			}
			if (fp)
			{
				while (!feof(fp) && fgets(name, 100, fp))
				{
					for (char* tmp = name; *tmp; tmp++)
					{
						if (*tmp == 10 || *tmp == 13)
						{
							*tmp = 0;
							break;
						}
					}
					std::string cmd(name);
					ProcessCommand(cmd);
				}
				fclose(fp);
				Status("Execution complete.");
			}
			else
			{
				Status("Could not read file.");
			}
		}
		else
		{
			Status("Could not read file.");
		}
	}
	else if (cmd.substr(0,8) == "football"){
		double l = 5.5f;
		double w = 3.5f;
		sscanf(cmd.c_str(), "football %lf %lf", &l, &w);
		this->football(l, w);
		Status("Football created.");
	}
	else if (cmd.substr(0, 4) == "move"){
		if (!mShapeStack.empty()){
			double x = 0;
			double y = 0;
			double z = 0;
			sscanf(cmd.c_str(), " move %lf %lf %lf", &x, &y, &z);
			this->move(x, y, z);
			Status("Shape moved.");
		}else{
			Status("Shape stack is empty.");
		}
	}
	else if (cmd.substr(0, 5) == "point"){
		double x = 0;
		double y = 0;
		double z = 0;
		sscanf(cmd.c_str(), " point %lf %lf %lf", &x, &y, &z);
		this->point(x, y, z);
		Status("Point created.");
	}
	else if (cmd == "pop"){
		if (!mShapeStack.empty()){
			this->pop();
			Status("Shape popped from stack.");
		}else{
			Status("Shape stack is empty.");
		}
	}
	else if (cmd.substr(0,4) == "rect")
	{
		double x = 5;
		double y = 5;
		sscanf(cmd.c_str(), " rect %lf %lf", &x, &y);
		this->rect(x, y);
		Status("Rectangle created.");
	}
	else if (cmd.substr(0, 3) == "rot"){
		if (!mShapeStack.empty()){
			double x = 0;
			double y = 0;
			double z = 0;
			sscanf(cmd.c_str(), " rot %lf %lf %lf", &x, &y, &z);
			this->rot(x, y, z);
			Status("Shape rotated.");
		}else{
			Status("Shape stack is empty.");
		}
	}
	else if (cmd.substr(0,6) == "saucer"){
		double r = 3;
		double t = 1;
		sscanf(cmd.c_str(), " saucer %lf %lf", &r, &t);
		this->saucer(r, t);
		Status("Saucer created.");
	}
	else if (cmd.substr(0,7) == "segment"){
		double x = 5;
		sscanf(cmd.c_str(), " segment %lf", &x);
		this->segment(x);
		Status("Segment created.");
	}
	else if (cmd.substr(0,6) == "sphere"){
		double rx = 5;
		double ry = 5;
		double rz = 5;
		int count = sscanf(cmd.c_str(), "sphere %lf %lf %lf", &rx, &ry, &rz);
		if (count == 1){
			this->sphere(rx);
		}else{
			this->sphere(rx, ry, rz);
		}
		Status("Sphere created.");
	}
	else if (cmd == "swap"){
		if (mShapeStack.size() >= 2){
			this->swap();
			Status("Shapes swapped.");
		}else{
			Status("Shape stack is empty.");
		}
	}
	else if (cmd == "sweep"){
		if (mShapeStack.size() >= 2){
			this->sweep();
			Status("Minkowski sum (sweep) created.");
		}else{
			Status("Minkowski sum requires two shapes.");
		}
	}
	else if (cmd == "wrap"){
		if (mShapeStack.size() >= 2){
			this->wrap();
			Status("Minkowski max (shrink wrap) created.");
		}else{
			Status("Minkowski max requires two shapes.");
		}
	}
	else{
		Status("Command not recognized.  Type 'help' for list of commands.");
	}
}

//////////////////////////////////////////////////////////////

void BodyBuilder::Status(const char* str)
{
	mStatusLine = str;
}

//////////////////////////////////////////////////////////////
