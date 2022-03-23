#include "./Observer.h"
// Здесь заинклюдим рандомчик и будем генерить радиyс и корды
Observer::Observer(){ //конструктор

}
void Observer::set_x(double new_x){
	coordinates[0] = new_x = x;
	return;
} 
void Observer::set_y(double new_y){
	coordinates[1] = new_y = y;
	return;
}
void Observer::set_z(double new_z){
	coordinates[2] = new_z = z;
	return;
}

void Observer::set_radius(double new_radius){ //радиус до центра фигуры
	radius = new_radius;
	return;
}

void Observer::replace(double* new_coordinates){ // Поменять сразу весь массив
	for(int i=0; i < 3; i++)
		coordinates[i] = new_coordinates[i];
	return;
} 
double Observer::get_x(){
	return x;
}
double Observer::get_y(){
	return y;
}
double Observer::get_z(){
	return z;
}
double* Observer::get_cords(){
	return coordinates;
}
double Observer::get_radius(){
	return radius;
}
double Observer::shoot(){ //Измерить расстояние
	//???????????????????????????????????
}