/*
 * Ellipse.cpp
 *
 *  Created on: 10.08.2017
 *      Author: ciesla
 */

#include "Ellipse.h"

double Ellipse::longSemiAxis;
double Ellipse::shortSemiAxis;
double Ellipse::neighbourListCellSize;
double Ellipse::voxelSize;

void Ellipse::calculateU(){
	this->u[0]  = cos(this->orientation); this->u[1]  = -sin(this->orientation);
	this->uT[0] = sin(this->orientation); this->uT[1] =  cos(this->orientation);
}

Ellipse::Ellipse() : AnisotropicShape2D(){
	this->a = Ellipse::longSemiAxis;
	this->b = Ellipse::shortSemiAxis;
	this->orientation = 0.0;
	this->calculateU();
}

Ellipse::~Ellipse() {
}

Ellipse & Ellipse::operator = (const Ellipse & el){
	this->a = el.a;
	this->b = el.b;
	this->orientation = el.orientation;
	for(unsigned char i=0; i<2; i++){
		this->u[i]  = el.u[i];
		this->uT[i] = el.uT[i];
	}
	return *this;
}

void Ellipse::initClass(const std::string &args){
	double ratio = std::stod(args);
	Ellipse::shortSemiAxis = sqrt(1.0/(M_PI*ratio));
	Ellipse::longSemiAxis = ratio*shortSemiAxis;
	Ellipse::neighbourListCellSize = 2*longSemiAxis;
	Ellipse::voxelSize = 1.4*shortSemiAxis;
}

Shape<2> * Ellipse::create(RND *rnd){
	Ellipse *el = new Ellipse();
	el->a = Ellipse::longSemiAxis;
	el->b = Ellipse::shortSemiAxis;
	el->orientation = rnd->nextValue()*2*M_PI;
	el->u[0]  = cos(el->orientation); el->u[1]  = sin(el->orientation);
	el->uT[0] = -sin(el->orientation); el->uT[1] = cos(el->orientation);
	return el;
}

double Ellipse::calculateF(double* r, double g){
	double d1 = (r[0]*this->u[0] + r[1]*this->u[1])  /this->a;
	double d2 = (r[0]*this->uT[0] + r[1]*this->uT[1]) /this->b;

	return 1 + g - d1*d1 - d2*d2;
}

void Ellipse::rotate(double* point, double alpha){
	double cosa = cos(alpha);
	double sina = sin(alpha);
	double x = point[0];
	double y = point[1];
	point[0] = x*cosa - y*sina;
	point[1] = x*sina + y*cosa;
}

int Ellipse::overlap(BoundaryConditions *bc, Shape *s) {
	Ellipse es = *((Ellipse *)s);
	double da[2];
	bc->getTranslation(da, this->position, es.position);
	es.translate(da);
	double d;
	d = (this->a/this->b - this->b/this->a)*sin(this->orientation - es.orientation);
	double g = 2 + d*d;
	double r[] = {this->position[0] - es.position[0], this->position[1] - es.position[1]};
	double f1 = this->calculateF(r, g);
	double f2 = es.calculateF(r, g);
	double psi = 4*(f1*f1-3*f2)*(f2*f2-3*f1) - (9-f1*f2)*(9-f1*f2);
	if (psi>0 && (f1<0 || f2<0)){
		return false;
	}
	return true;
}

double Ellipse::getVolume() {
	return M_PI*this->a*this->b;
}

int Ellipse::pointInside(BoundaryConditions *bc, double* da) {
	double ta[2];
	double tmp[2];

	tmp[0] = da[0]; tmp[1] = da[1];
	bc->getTranslation(ta, this->position, tmp);
	da[0] += ta[0];
	da[1] += ta[1];

	da[0] -= this->position[0];
	da[1] -= this->position[1];

	Ellipse::rotate(da, -this->orientation);

	double dx = da[0]/(this->a+this->b);
	double dy = da[1]/(2*this->b);
	return (dx*dx+dy*dy < 1);
}


int Ellipse::pointInside(BoundaryConditions *bc, double* da, double angleFrom, double angleTo) {


	// to be inside an exclusion zone the point have to bee inside all exclusion zones between angleFrom and angleTo

	Ellipse *ellTmp = new Ellipse();
	ellTmp->position[0] = da[0]; ellTmp->position[1] = da[1];

	// checking if the point is inside the exclusion zone for angleFrom
	ellTmp->orientation = angleFrom;
	ellTmp->calculateU();
	if( !this->overlap(bc, ellTmp) ){
		delete ellTmp;
		return false;
	}

	// checking if the point is inside the exclusion zone for angleTo
	ellTmp->orientation = angleTo;
	ellTmp->calculateU();
	if ( !this->overlap(bc, ellTmp)){
		delete ellTmp;
		return false;
	}

	// transforming coordinate system to one connected with "this" ellipse
	double point[2];
	point[0] = da[0] - this->position[0];
	point[1] = da[1] - this->position[1];
	angleFrom -= this->orientation;
	angleTo -= this->orientation;

	// TODO



	delete ellTmp;
	return true;

}



/*
	private void drawEllipse(Graphics g, double scale) {
		Polygon p = new Polygon();
		final int points = 200;

		double[] da = new double[2];
		for(int i=0; i<points; i++){
			double t = i*2*Math.PI/points;
			da[0] = a*Math.cos(t);
			da[1] = b*Math.sin(t);
			Calculator.rotate(da, this.orientation);
			da[0] += this.position[0];
			da[1] += this.position[1];
			p.addPoint((int)(da[0]*scale), (int)(da[1]*scale));
		}
		g.fillPolygon(p);

//		g.setColor(Color.BLACK);
//		byte[] no = String.valueOf(this.no).getBytes();
//		g.drawBytes(no, 0, no.length, (int)((position[0])*scale), (int)((position[1])*scale));
//		g.setColor(Color.RED);

	}

	private void drawExclusion(Graphics g, double scale) {
		Polygon p = new Polygon();
		final int points = 200;

		double[] da = new double[2];
		for(int i=0; i<points; i++){
			double t = i*2*Math.PI/points;
			da[0] = (a+b)*Math.cos(t);
			da[1] = 2*b*Math.sin(t);
			Calculator.rotate(da, this.orientation);
			da[0] += this.position[0];
			da[1] += this.position[1];
			p.addPoint((int)(da[0]*scale), (int)(da[1]*scale));
		}
		g.drawPolygon(p);

//		g.setColor(Color.BLACK);
//		byte[] no = String.valueOf(this.no).getBytes();
//		g.drawBytes(no, 0, no.length, (int)((position[0])*scale), (int)((position[1])*scale));
//		g.setColor(Color.RED);

	}

	public void draw(Graphics g, double scale) {
		this.drawExclusion(g, scale);
		this.drawEllipse(g, scale);
	}

	@Override
	public void shapeSpecificVoxelAnalysis(Surface surf, VoxelList voxels) {
		// nothing to do
	}

*/

double Ellipse::getNeighbourListCellSize() {
	return Ellipse::neighbourListCellSize;
}

double Ellipse::getVoxelSize() {
	return Ellipse::voxelSize;
}

std::string Ellipse::toPovray() const{
	//TODO
	return "";
}

void Ellipse::store(std::ostream &f) const{
	Shape::store(f);
	f.write((char *)(&this->a), sizeof(double));
	f.write((char *)(&this->b), sizeof(double));
	f.write((char *)(&this->orientation), sizeof(double));
}

void Ellipse::restore(std::istream &f){
	Shape::restore(f);
	f.read((char *)(&this->a), sizeof(double));
	f.read((char *)(&this->b), sizeof(double));
	f.read((char *)(&this->orientation), sizeof(double));
	this->calculateU();

}
