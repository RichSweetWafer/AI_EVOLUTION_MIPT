#include "Observer.h"
#include "./random_utilities.h"
#include <iostream>

//\/!!!!TEST RUN!!!\/\/


int main(){
	double* center =  new double [3] {0.0, 0.0, 0.0};
	Observer agent = Observer(center, 1000.0);
	std::cout << agent.get_x() << std::endl;
	std::cout << agent.get_y() << std::endl;
	std::cout << agent.get_z() << std::endl;
	return 0;
}

///!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!