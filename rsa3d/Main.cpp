#include "PackingGenerator.h"
#include "Shape.h"
#include <fstream>
#include <iostream>

void toPovRay(char *filename, double size, PackingGenerator &pg){
	std::fstream file(filename, std::ios::out);

	file << "#include \"colors.inc\"" << std::endl;
	file << "background { color White }" << std::endl;
	file << "camera { orthographic location <" << size/2 << ", " << size/2 << ", " << size << "> look_at  <" << size/2 << ", " << size/2 << ",  0> }" << std::endl;
	file << "light_source { < 1000.0, 1000.0, 1000.0> color White shadowless parallel point_at <" << size/2 << ", " << size/2 << ",  0>}" << std::endl;
	file << "#declare layer=union{" << std::endl;

	file << "  polygon {4, <0, 0, 0.0>, <0, " << size << ", 0.0>, <" << size << ", " << size << ", 0.0>, <" << size << ", 0, 0.0>  texture { finish { ambient 1 diffuse 0 } pigment { color Gray} } }" << std::endl;
	for(Shape *s : *pg.getPacking()){
		double *da = s->getPosition();
		file << "  sphere { <" << da[0] << ", " << da[1] << ", 0.0>, " << 0.56418958354775628 << std::endl;
		file << "    texture { pigment { color Red } }" << std::endl;
		file << "  }" << std::endl;
		std::cout << s->no << "\t";
	}


	file << "}" << std::endl;
	file << "#declare result=union{" << std::endl;
	file << "  object { layer }" << std::endl;
	file << "}" << std::endl;
	file << "object{ result	rotate x*360*clock }" << std::endl;

	file.close();
}


int main(int argc, char **argv){
	Parameters params;
	PackingGenerator pg(1, &params);
	pg.run();

	toPovRay("surface.pov", 100.0, pg);
/*
	for(Shape *s : *pg.getPacking()){
		double *da = s->getPosition();
		std::cout << s->no << "\t";
		for(int i=0; i<2; i++)
			std::cout << da[i] << "\t";
		std::cout << std::endl;
	}
*/
}


