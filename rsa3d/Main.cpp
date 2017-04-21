#include "PackingGenerator.h"
#include "Shape.h"
#include <iostream>

int main(int argc, char **argv){
	Parameters params;
	PackingGenerator pg(1, &params);
	pg.run();

	for(Shape *s : *pg.getPacking()){
		double *da = s->getPosition();
		std::cout << s->no << "\t";
		for(int i=0; i<2; i++)
			std::cout << da[i] << "\t";
		std::cout << std::endl;
	}
}

