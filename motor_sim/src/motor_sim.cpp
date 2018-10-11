//============================================================================
// Name        : prova2.cpp
// Author      :
// Version     :
// Copyright   : Your copyright notice
// Description : Hello World in C++, Ansi-style
//============================================================================

#include <Windows.h>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <iostream>

using namespace std;

float Ralphae;
float Romegae;
float Rthetae;
float Salphae;
float Somegae_prefilter;
float Somegae_postfilter;
float Somegae;
float Somegae_windmill;
float Somegae_tod;
float Somegae_deltalim;
float Sthetae;
float Samplie;
float Samplie_clamped;
float Samplie_alt;
float Lalphae;
float Walphae;
float tod_f;
float tod;
struct tod_sampling_t {
  float value;
  float time;
  float interval;
  float factor;
  float tod_timeconstant;
  float outfilter_timeconstant;
} tod_sampling;
float t;
float t_inc;
int indexx;
int gve;
int gvi;
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
struct windmill_st {
  float integral;
  float proportional;
  float out;
  float KI;
  float KP;
} windmill;
float a, b, c, max_triangle;  // per clamp tensione
float vbatt;
struct delta_lim_st {
  float delta_max;
  float delta_lim;
  float KI;
  float integral;
} delta_lim;
float pow_windmill;
FILE* file1;

const float PI = 3.14159F;
const float TUP = 16;
const float TDOWN = 12;
const float TSTILL = 10;
const float EXTENDED_TIME = 15;
const float T_WINDMILL_UP = 20;
const float T_WINDMILL_STILL = 5;
const float T_WINDMILL_DOWN = 20;
const float T_WINDMILL_AFTER = 20;
const float DT = 50e-6;
const float FLUX = 1e-2;
const float R = 33.8e-3;
const float L = 223e-6;
const float PP = 4;
const float JE = (5.46e-3 / PP / PP);
const float ROMEGAE_MAX = (250.0 * 2 * PI);
const float ROMEGAE_TH1 = (60.0 * 2 * PI);
const float SOMEGAE_MAX = (550.0 * 2 * PI);
const float PMAX = 550.0;
const float PDEC = 200.0;
const float SOMEGAE_TIMECONST = 1.6;
const float SAMPLIE_OFFSET = 0.61;

#define DEGREE_TO_RADIANS(x) (x / 180.0 * PI)

int main() {
  Ralphae = 0;
  Romegae = 0;
  Rthetae = 0;
  Salphae = 0;
  Somegae = 0;
  Sthetae = 0;
  tod_f = 0;
  tod = 0;
  ix = 0;
  iy = 0;
  gamma.iy_last = 0;
  t_inc = 1 / Somegae / 100;
  gve = 0;
  gvi = 0;
  tod_sampling.time = 0;
  tod_sampling.value = 0;
  tod_sampling.interval = 1e-3;
  tod_sampling.factor = 10;
  tod_sampling.tod_timeconstant = 1;
  tod_sampling.outfilter_timeconstant = 0.1;
  gamma.integral = 0;
  gamma.KI = 400;

  i_min.THRESHOLD = 6;
  i_min.integral = 0;
  i_min.KI = 1;

  vbatt = 24.0;

  delta_lim.KI = 4000;
  delta_lim.integral = 0;

  windmill.KI = 500;
  windmill.KP = 10;
  windmill.integral = 0;
  windmill.out = 0;

  pow_windmill = 0;

  indexx = 0;

  file1 = fopen("pippo.txt", "w");
  if (file1 == NULL) {
    return -1;
  }

  for (t = 0; t < TUP + TSTILL + TDOWN + EXTENDED_TIME + T_WINDMILL_UP +
                      T_WINDMILL_STILL + T_WINDMILL_DOWN + T_WINDMILL_AFTER;
       t += t_inc) {
    // Calcolo incremento temporale.
    // Massimo 1/100 radiante o 10 ms
    if (Somegae != 0) t_inc = 1 / Somegae / 100;
    if (t_inc > 0.01) t_inc = 0.01;

    // Break simulation if control is broken.
    // Sync is lost
    float q1 = Sthetae - Rthetae + 2 * PI * gve - PI / 2;
    if (q1 > 2 * PI || q1 < -2 * PI) {
      break;
    }
    // Control is windmilling and voltage is at max
    if (gamma.phi > DEGREE_TO_RADIANS(75.0)) {
      if (Samplie * sqrt(3) > vbatt * 700.0 / 512) {
        break;
      }
    }

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
    if (t < TUP) {
      Salphae = ROMEGAE_MAX / TUP;
    } else if (t > (TUP + TSTILL) && (t < (TUP + TSTILL + TDOWN))) {
      Salphae = (ROMEGAE_TH1 + 10 - ROMEGAE_MAX) / TDOWN;
    } else {
      Salphae = 0;
    }

    // potenza descrescente per simulare windmill
    if (t > TUP + TDOWN + TSTILL + EXTENDED_TIME &&
        t < TUP + TDOWN + TSTILL + EXTENDED_TIME + T_WINDMILL_UP) {
      pow_windmill += (PDEC / T_WINDMILL_UP) * t_inc;
    }
    if (t > (TUP + TDOWN + TSTILL + EXTENDED_TIME + T_WINDMILL_UP +
             T_WINDMILL_STILL) &&
        t < (TUP + TDOWN + TSTILL + EXTENDED_TIME + T_WINDMILL_UP +
             T_WINDMILL_STILL + T_WINDMILL_DOWN)) {
      if (pow_windmill > 0) {
        pow_windmill -= (PDEC / T_WINDMILL_DOWN) * t_inc;
      } else {
        pow_windmill = 0;
      }
    }

    Somegae_prefilter += Salphae * t_inc;
    // Clamp velocit'a max statore
    if (Somegae_prefilter > SOMEGAE_MAX) {
      Somegae_prefilter = SOMEGAE_MAX;
    }
    // filtro passa basso su velocita' statore (sigmoide)
    Somegae_postfilter +=
        (Somegae_prefilter - Somegae_postfilter) * t_inc / SOMEGAE_TIMECONST;
    // Applicazione windmill
    Somegae = Somegae_postfilter + windmill.out;
    // Correzione delta lim
    Somegae_deltalim = Somegae - delta_lim.integral;
    // Applicazione tod
    Somegae_tod = Somegae_deltalim - tod_sampling.value;
    Sthetae += Somegae_tod * t_inc;
    // Rotore
    Romegae += Ralphae * t_inc;
    Rthetae += Romegae * t_inc;

    // Troncamento angoli e calcolo di overflow
    if (Sthetae > 2 * PI) {
      Sthetae -= 2 * PI;
      gve++;
      gvi++;
    }
    if (Rthetae > 2 * PI) {
      Rthetae -= 2 * PI;
      gve--;
    }

    // Ampiezza tensione statorica
    Samplie_alt = 1;
    Samplie_alt +=
        sqrt(pow(Somegae * L * i, 2) + pow((R * i + Somegae * FLUX), 2));

    Samplie = SAMPLIE_OFFSET + Somegae * FLUX + gamma.integral;

    // Overmodulation clamp (24V battery)
    a = sin(Sthetae);
    b = sin(Sthetae + 2.0 / 3 * PI);
    c = sin(Sthetae - 2.0 / 3 * PI);
    max_triangle = 0;
    if (a - b > max_triangle) max_triangle = a - b;
    if (b - a > max_triangle) max_triangle = b - a;
    if (a - c > max_triangle) max_triangle = a - c;
    if (c - a > max_triangle) max_triangle = c - a;
    if (b - c > max_triangle) max_triangle = b - c;
    if (c - b > max_triangle) max_triangle = c - b;
    if (max_triangle * Samplie > vbatt) {
      Samplie_clamped = vbatt / max_triangle;
      // gamma.integral = Samplie - (SAMPLIE_OFFSET + Somegae_tod * FLUX);
    } else {
      Samplie_clamped = Samplie;
    }

    // Aggiornamento corrente
    dix = (Samplie_clamped * cos(Sthetae) -
           Romegae * FLUX * cos(Rthetae + PI / 2) - R * ix) /
          L * t_inc;
    diy = (Samplie_clamped * sin(Sthetae) -
           Romegae * FLUX * sin(Rthetae + PI / 2) - R * iy) /
          L * t_inc;
    ix += dix;
    iy += diy;
    i = sqrt(ix * ix + iy * iy);

    // Decelerazione dovuta al carico
    Lalphae = PMAX * pow(Romegae, 2) / pow(ROMEGAE_MAX, 3) / JE;
    // Accelerazione dovuta a windmill
    if (Romegae > 0) {
      Walphae = pow_windmill / Romegae / JE;
    }
    // Accelerazione dovuta a corrente statorica
    Ralphae = 3 / 2 * FLUX * (cos(Rthetae) * iy - sin(Rthetae) * ix) / JE;

    Ralphae = Ralphae - Lalphae + Walphae;

    // TOD
    tod = (Sthetae - Rthetae + 2 * PI * gve) - tod_f;
    tod_f += tod * t_inc / tod_sampling.tod_timeconstant;
    if (t - tod_sampling.time > tod_sampling.interval) {
      tod_sampling.time = t;
      if (Romegae > ROMEGAE_TH1) {
        tod_sampling.value = tod * tod_sampling.factor;
      } else {
        tod_sampling.value = 0;
        // this will make a smooth tod kick in
        tod_f = Sthetae - Rthetae + 2 * PI * gve;
      }
    }

    // Corrente minima
    if (i < i_min.THRESHOLD || i_min.integral > 0) {
      i_min.integral += (i_min.THRESHOLD - i) * t_inc * i_min.KI;
    }

    // gamma
    // detect zero cross
    if (iy >= 0 && gamma.iy_last < 0) {
      gvi--;  // overflow v-i difference
      gamma.phi = Sthetae + 2 * PI * gvi;
      gamma.delta = Somegae * L * i / Samplie_clamped;
      gamma.gamma = gamma.delta - gamma.phi;
      // gamma.gamma = -sin(Rthetae + PI / 2);
      if (Romegae > ROMEGAE_TH1) {
        gamma.integral += gamma.gamma * t_inc * gamma.KI + i_min.integral;
      } else {
        gamma.integral = 0;
      }
    }
    gamma.iy_last = iy;

    // Delta Angle max value management
    delta_lim.delta_max = atan(Somegae * L / R);
    delta_lim.delta_lim = 2.0 / 3 * delta_lim.delta_max;
    if (Romegae > ROMEGAE_TH1) {
      float q = Sthetae - Rthetae + 2 * PI * gve - PI / 2 - delta_lim.delta_lim;
      if (q > 0 || delta_lim.integral > 0) {
        delta_lim.integral += q * t_inc * delta_lim.KI;
      }
    }

    // Windmill
    float q = gamma.phi - DEGREE_TO_RADIANS(75);
    bool c1 = 0 < q && ROMEGAE_TH1 < Romegae;
    bool c2 = 0 < windmill.out && q < 0;
    if (c1 || c2) {
      windmill.integral += q * t_inc * windmill.KI;
      windmill.proportional = q * windmill.KP;
      if (windmill.proportional < 0) {
        windmill.proportional = 0;
      }
      windmill.out = windmill.integral + windmill.proportional;
    }
    else {
    	windmill.out = 0;
    	windmill.integral = 0;
    }

    // Stampa variabili
    if (indexx++ % 1000 == 0) {
      char sf[2048] = "";
      char s1[2048] = "";
      if (indexx == 1) {
        fprintf(file1,
                "t\t"
                "Samplie\t"
                "Samplie_alt\t"
                "Lalphae\t"
                "Ralphae\t"
                "Romegae\t"
                "Rthetae\t"
                "Salphae\t"
                "Somegae_tod\t"
                "(Sthetae-Rthetae+2*PI*gve-PI/2)\t"
                "ix\t"
                "iy\t"
                "i\t"
                "tod\t"
                "delta_lim.delta_lim\t"
                "delta_lim.delta_max\t"
                "pow_windmill\t"
                "gamma.phi\t"
                "windmill.out\t"
                "delta_lim.integral\t"
                "gamma.integral\t"
                "i_min.integral\t"
                "windmill.out\t"
                "Romegae*FLUX\t"
                "Somegae\t"
                "gvi\t");
      }
      sprintf(s1, "%.2f\t ", /*01*/ t);
      strcat(sf, s1);
      sprintf(s1, "%.2f\t ", /*02*/ Samplie);
      strcat(sf, s1);
      sprintf(s1, "%.2f\t ", /*03*/ Samplie_alt);
      strcat(sf, s1);
      sprintf(s1, "%.2f\t ", /*04*/ Lalphae);
      strcat(sf, s1);
      sprintf(s1, "%.2f\t ", /*05*/ Ralphae);
      strcat(sf, s1);
      sprintf(s1, "%.2f\t ", /*06*/ Romegae);
      strcat(sf, s1);
      sprintf(s1, "%.2f\t ", /*07*/ Rthetae);
      strcat(sf, s1);
      sprintf(s1, "%.2f\t ", /*08*/ Salphae);
      strcat(sf, s1);
      sprintf(s1, "%.2f\t ", /*09*/ Somegae_tod);
      strcat(sf, s1);
      sprintf(s1, "%.2f\t ",
              /*10*/ (Sthetae - Rthetae + 2 * PI * gve - PI / 2));
      strcat(sf, s1);
      sprintf(s1, "%.2f\t ", /*11*/ ix);
      strcat(sf, s1);
      sprintf(s1, "%.2f\t ", /*12*/ iy);
      strcat(sf, s1);
      sprintf(s1, "%.2f\t ", /*13*/ i);
      strcat(sf, s1);
      sprintf(s1, "%.2f\t ", /*14*/ tod);
      strcat(sf, s1);
      sprintf(s1, "%.2f\t ", /*15*/ delta_lim.delta_lim);
      strcat(sf, s1);
      sprintf(s1, "%.2f\t ", /*16*/ delta_lim.delta_max);
      strcat(sf, s1);
      sprintf(s1, "%.2f\t ", /*17*/ pow_windmill);
      strcat(sf, s1);
      sprintf(s1, "%.2f\t ", /*18*/ gamma.phi);
      strcat(sf, s1);
      sprintf(s1, "%.2f\t ", /*19*/ windmill.out);
      strcat(sf, s1);
      sprintf(s1, "%.2f\t ", /*20*/ delta_lim.integral);
      strcat(sf, s1);
      sprintf(s1, "%.2f\t ", /*21*/ gamma.integral);
      strcat(sf, s1);
      sprintf(s1, "%.2f\t ", /*22*/ i_min.integral);
      strcat(sf, s1);
      sprintf(s1, "%.2f\t ", /*23*/ windmill.out);
      strcat(sf, s1);
      sprintf(s1, "%.2f\t ", /*24*/ Romegae * FLUX);
      strcat(sf, s1);
      sprintf(s1, "%.2f\t ", /*25*/ Somegae);
      strcat(sf, s1);
      sprintf(s1, "%d\t ", /*26*/ gvi);
      strcat(sf, s1);
      strcat(sf, "\n");
      fprintf(file1, sf);
      // Sleep(50);
      // fflush(stdout);
    }
  }
  fclose(file1);
  return 0;
}
