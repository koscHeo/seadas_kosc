void l_sun(int iyr, int iday, double sec, float sunr[3], float *rs);
void sun2000(int iyr, int iday, double sec, float sun[3], float *rs);
void gha2000(int iyr, double day, double *gha);
void ephparms(double t, double *xls, double *gs, double *xlm, double *omega);
void nutate(double t, double xls, double gs, double xlm, double omega,
            double *dpsi, double *eps);
long jd(int year, int month, int day);
