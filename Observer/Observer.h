//Наблюдатель появляется в круге радиуса N от центра фигуры.

class Observer{
public:
	Observer(double* center, double radius); //конструктор
	~Observer(){
		delete coordinates;
	}
	void set_x(double new_x); 
	void set_y(double new_y);
	void set_z(double new_z);
	void replace(double* new_coordinates); // Поменять сразу весь массив 
	double get_x();
	double get_y();
	double get_z();
	double* get_cords();
	double get_radius(); //Насколько это интересно наблюдателю? Может отдельный класс среды?
	double shoot(); //Измерить расстояние

private:
	double* coordinates;
	double radius;
	double* center;
	double x;
	double y;
	double z;
};