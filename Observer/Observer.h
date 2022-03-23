//Наблюдатель появляется в круге радиуса N от центра фигуры.

class Observer{
public:
	Observer(); //конструктор
	void set_x(double new_x); 
	void set_y(double new_y);
	void set_z(double new_z);
	void set_radius(double new_radius); //радиус до центра фигуры
	void replace(double* new_coordinates); // Поменять сразу весь массив 
	double get_x();
	double get_y();
	double get_z();
	double* get_cords();
	double get_radius();
	double shoot(); //Измерить расстояние

private:
	double coordinates[3];
	double radius;
	//double center; Как он поймет где фигура?
	double x;
	double y;
	double z;
};