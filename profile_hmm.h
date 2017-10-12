
typedef struct{
  int *flagX;
  int *flagY;
  int Nx;
  int Ny;
  float (*scoreX)[N_PEPT];
  float (*scoreY)[N_PEPT];
}profile;

void init_prof_HMM(int len);
void calc_prof_posterior(int len);
