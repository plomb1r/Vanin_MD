// Кооординаты (массивы от числа частиц, где i - номер частицы)
double* coordx;
double* coordy;
double* coordz;

// Скорости (массивы от числа частиц, где i - номер частицы)
double* vx;	 // скорости 
double* vy;
double* vz;


double U(double r); // потенциальная энергия(скаляр)


// Силы взаимодействия (массивы от числа частиц и вспомогательный скаляр)
double* Fx;
double* Fy;
double* Fz;

double F(double r);