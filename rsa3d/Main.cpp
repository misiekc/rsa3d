#include "PackingGenerator.h"
#include "Shape.h"
#include "shapes/Sphere.h"
#include <fstream>
#include <iostream>

void toPovRay(std::string filename, double size, std::vector<Shape *> *packing){
	std::ofstream file(filename);

	file << "#include \"colors.inc\"" << std::endl;
	file << "background { color White }" << std::endl;
	file << "camera { orthographic location <" << size/2 << ", " << size/2 << ", " << (2*size) << "> look_at  <" << size/2 << ", " << size/2 << ",  0> }" << std::endl;
	file << "light_source { < 1000.0, 1000.0, 1000.0> color White shadowless parallel point_at <" << size/2 << ", " << size/2 << ",  0>}" << std::endl;
	file << "#declare layer=union{" << std::endl;

//	file << "  polygon {4, <0, 0, 0.0>, <0, " << size << ", 0.0>, <" << size << ", " << size << ", 0.0>, <" << size << ", 0, 0.0>  texture { finish { ambient 1 diffuse 0 } pigment { color Gray} } }" << std::endl;
	for(Shape *s : *packing){
		file << s->toPovray();
	}


	file << "}" << std::endl;
	file << "#declare result=union{" << std::endl;
	file << "  object { layer }" << std::endl;
	file << "}" << std::endl;
	file << "object{ result	rotate x*360*clock }" << std::endl;

	file.close();
}

void toFile(std::string filename, std::vector<Shape *> *packing){
	std::ofstream file(filename, std::ios::binary);
	for(Shape *s : *packing){
		s->store(file);
	}
	file.close();
}

std::vector<Shape *> * fromFile(std::string filename){
	std::ifstream file(filename, std::ios::binary);
	std::vector<Shape *> * v = new std::vector<Shape *>;

	while(!file.eof()){
		Shape *s = Sphere::restore(file);
		v->push_back(s);
	}
	v->pop_back();

	file.close();
	return v;
}


int main(int argc, char **argv){
	Parameters params(argv[1]);
	PackingGenerator pg(1, &params);
	pg.run();

	toPovRay("cubs.pov", 10.0, pg.getPacking());
//	toFile("surf.bin", pg.getPacking());
}

/*
int main(int argc, char **argv){

	std::vector<Shape *> *packing = fromFile("surf.bin");
	toPovRay("surf2.pov", 100.0, packing);
}
*/
