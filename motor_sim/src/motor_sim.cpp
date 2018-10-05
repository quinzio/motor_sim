//============================================================================
// Name        : prova2.cpp
// Author      : 
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <iostream>
#include <cmath>
#include <Windows.h>
#include <cstdio>
#include <cstdlib>

using namespace std;

float Ralphae;
float Romegae;
float Rthetae;
float Salphae;
float Somegae_prefilter;
float Somegae;
float Somegae_tod;
float Sthetae;
float Samplie;
float Samplie_alt;
float Lalphae;
float tod_f;
float tod;
struct tod_sampling_t {
	float value;
	float time;
	float interval;
	float factor;
	float tod_timeconstant;
} tod_sampling;
float t;
float t_inc;
int indexx;
int gd;
float ix;
float iy;
struct gamma_t {
	float iy_last;
	float phi;
	float gamma;
	float delta;
	float integral;
	float KI;
} gamma;
float dix;
float diy;
float z;
float i;
struct i_min_st {
	float THRESHOLD;
	float integral;
	float KI;
} i_min;
FILE * file1;

#define PI 3.14159F
#define TTOT 16
#define EXTENDED_TIME 20
#define DT 50e-6
#define FLUX 1e-2
#define R 33.8e-3
#define L 223e-6
#define PP 4
#define JE (5.46e-3 / PP / PP)
#define ROMEGAE_MAX (250.0 * 2 * PI)
#define ROMEGAE_TH1 (60.0 * 2 * PI)
#define SOMEGAE_MAX (66.0 * 2 * PI)
#define PMAX 500.0
#define SOMEGAE_TIMECONST 0.6

int main() {

	Ralphae = 0;
	Romegae = 0;
	Rthetae = 0;
	Salphae = 0;
	Somegae = 0;
	Sthetae = PI / 2;
	tod_f = 0;
	tod = 0;
	ix = 0;
	iy = 0;
	gamma.iy_last = 0;
	t_inc = 1 / Somegae / 100;
	gd = 0;
	tod_sampling.time = 0;
	tod_sampling.value = 0;
	tod_sampling.interval = 1e-3;
	tod_sampling.factor = 10;
	tod_sampling.tod_timeconstant = 1;
	gamma.integral = 0;
	gamma.KI = 400;

	i_min.THRESHOLD = 7;
	i_min.integral = 0;
	i_min.KI = 2;

	indexx = 0;

	file1 = fopen("pippo.txt", "w");
	if (file1 == NULL) {
		return -1;
	}

	for (t = 0; t < TTOT + EXTENDED_TIME; t += t_inc) {

		// Calcolo incremento temporale.
		// Massimo 1/100 radiante o 10 ms
		if (Somegae != 0)
			t_inc = 1 / Somegae / 100;
		if (t_inc > 0.01)
			t_inc = 0.01;

		// Eq. alle differenze moto rotore e statore
		// Accelerazione statore "quadratica"
		/*
		 if (t < TTOT) {
		 Salphae = (2.0 * ROMEGAE_MAX / TTOT * (1 - t / TTOT));
		 } else {
		 Salphae = 0;
		 }
		 */
		// Accelerazione statore "lineare"
		if (t < TTOT) {
			Salphae = ROMEGAE_MAX / TTOT;
		} else {
			Salphae = 0;
		}
		Somegae_prefilter += Salphae * t_inc;
		// Clamp velocit'a max statore
		if (Somegae_prefilter > SOMEGAE_MAX) {
			Somegae_prefilter = SOMEGAE_MAX;
		}
		// filtro passa basso su velocita' statore (sigmoide)
		Somegae += (Somegae_prefilter - Somegae) * t_inc / SOMEGAE_TIMECONST;
		// Applicazione tod
		Somegae_tod = Somegae - tod_sampling.value;
		Sthetae += Somegae_tod * t_inc;
		// Rotore
		Romegae += Ralphae * t_inc;
		Rthetae += Romegae * t_inc;

		// Troncamento angoli e calcolo di overflow
		if (Sthetae > 2 * PI) {
			Sthetae -= 2 * PI;
			gd++;
		}
		if (Rthetae > 2 * PI) {
			Rthetae -= 2 * PI;
			gd--;
		}

		// Ampiezza tensione statorica
		Samplie_alt = 1;
		Samplie_alt += sqrt(
				pow(Somegae * L * i, 2) + pow(( R * i + Somegae * FLUX), 2));

		Samplie = 3.1 + Somegae_tod * FLUX + gamma.integral;

		// Aggiornamento corrente
		dix = (Samplie * cos(Sthetae) - Romegae * FLUX * cos(Rthetae + PI / 2)
				- R * ix) / L * t_inc;
		diy = (Samplie * sin(Sthetae) - Romegae * FLUX * sin(Rthetae + PI / 2)
				- R * iy) / L * t_inc;
		ix += dix;
		iy += diy;
		i = sqrt(ix * ix + iy * iy);

		// Decelerazione dovuta al carico
		Lalphae = PMAX * pow(Romegae, 2) / pow(ROMEGAE_MAX, 3) / JE;
		// Accelerazione dovuta a corrente statorica - decelerazione
		Ralphae = 3 / 2 * FLUX * (cos(Rthetae) * iy - sin(Rthetae) * ix) / JE
				- Lalphae;

		// TOD
		tod = (Sthetae - Rthetae + 2 * PI * gd) - tod_f;
		tod_f += tod * t_inc / tod_sampling.tod_timeconstant;
		if (t - tod_sampling.time > tod_sampling.interval) {
			tod_sampling.time = t;
			if (Romegae > ROMEGAE_TH1) {
				tod_sampling.value = tod * tod_sampling.factor;
			} else {
				tod_sampling.value = 0;
			}
		}

		// Corrente minima
		if (i < i_min.THRESHOLD || i_min.integral > 0) {
			i_min.integral += (i_min.THRESHOLD - i) * t_inc * i_min.KI;
		}

		// gamma
		// detect zero cross
		if (iy >= 0 && gamma.iy_last <= 0) {
			gamma.phi = Sthetae;
			gamma.delta = Somegae * L * i / Samplie;
			gamma.gamma = gamma.delta - gamma.phi;
			//gamma.gamma = -sin(Rthetae + PI / 2);
			if (Romegae > ROMEGAE_TH1) {
				gamma.integral += gamma.gamma * DT * gamma.KI + i_min.integral;
			} else {
				gamma.integral = 0;
			}

		}
		gamma.iy_last = iy;

		// Stampa variabili
		if (indexx++ % 100 == 0) {
			fprintf(file1,
					"%.2f\t %.2f\t %.2f\t %.2f\t %.2f\t %.2f\t %.2f\t %.2f\t %.2f\t %.2f\t %.2f\t %.2f\t %.2f\t %.5f\t %.2f\t \n",
					/*1*/t,
					/*2*/Samplie,
					/*3*/Samplie_alt,
					/*4*/Lalphae,
					/*5*/Ralphae,
					/*6*/Romegae,
					/*7*/Rthetae,
					/*8*/Salphae,
					/*9*/Somegae_tod,
					/*10*/(Sthetae - Rthetae + 2 * PI * gd),
					/*11*/ix,
					/*12*/iy,
					/*13*/i,
					/*14*/tod,
					/*15*/gamma.gamma);
			//Sleep(50);
			//fflush(stdout);
		}
	}
	fclose(file1);
	return 0;
}
