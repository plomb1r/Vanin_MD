#pragma once
void start_cond_two_particles() {
	const double r_vec_1[] = { 0.75, 0.75, 0.5 };
	const double r_vec_2[] = { 1.25, 0.75, 0.5 };
	const double v_vec_1[] = { 1.0, 1.0, 0.0 };
	const double v_vec_2[] = { -1.0, 1.0, 0.0 };
#pragma region Инициализация НУ для координат
	coordx[0] = r_vec_1[0];
	coordy[0] = r_vec_1[1];
	coordz[0] = r_vec_1[2];
	coordx[1] = r_vec_2[0];
	coordy[1] = r_vec_2[1];
	coordz[1] = r_vec_2[2];
#pragma endregion
#pragma region Инициализация НУ для скоростей
	vx[0] = v_vec_1[0];
	vy[0] = v_vec_1[1];
	vz[0] = v_vec_1[2];
	vx[1] = v_vec_2[0];
	vy[1] = v_vec_2[1];
	vz[1] = v_vec_2[2];
#pragma endregion
};
