
/////////////// Output Product Arrays ///////////////
float faraday_deg[MAXCYC*NUMBER_OF_BEAMS];
float sss_aquarius_nolc[MAXCYC*NUMBER_OF_BEAMS];
float sss_aquarius[MAXCYC*NUMBER_OF_BEAMS];
float sss_aquarius_adj[MAXCYC*NUMBER_OF_BEAMS];
float sss_aquarius_err_ran[MAXCYC*NUMBER_OF_BEAMS];
float sss_aquarius_err_sys[MAXCYC*NUMBER_OF_BEAMS];
float vol_sm[MAXCYC*NUMBER_OF_BEAMS];
//float d_sss_aquarius_land_corr[MAXCYC*NUMBER_OF_BEAMS];
float ta_earth[NUMBER_OF_BEAMS*3];
float tb_toi[NUMBER_OF_BEAMS*3];
float tb_toa_nolc[NUMBER_OF_BEAMS*2];
float tb_toa[NUMBER_OF_BEAMS*2];
float tb_sur[MAXCYC*NUMBER_OF_BEAMS*2];
float tb_sur_nolc[MAXCYC*NUMBER_OF_BEAMS*2];
float tb_sur_rc_nolc[MAXCYC*NUMBER_OF_BEAMS*2];
float tb_sur_rc[MAXCYC*NUMBER_OF_BEAMS*2];
float ta_expected[MAXCYC*NUMBER_OF_BEAMS*3];
float ta_expected_hhh[MAXCYC*NUMBER_OF_BEAMS*3];
float tb_expected[MAXCYC*NUMBER_OF_BEAMS*2];
float tb_expected0[MAXCYC*NUMBER_OF_BEAMS*2];
float tb_con_nolc[MAXCYC*NUMBER_OF_BEAMS];
float tb_con[MAXCYC*NUMBER_OF_BEAMS];
float tagal_ref_GO[NUMBER_OF_BEAMS*3];
float dtagal_ref[NUMBER_OF_BEAMS*3];

float density[MAXCYC*NUMBER_OF_BEAMS];
float spiciness[MAXCYC*NUMBER_OF_BEAMS];

float delta_tb_sst_wspd[MAXCYC*NUMBER_OF_BEAMS*2];

// Scatterometer output arrays
float scat_ant[NUMBER_OF_BEAMS*4];
float scat_toa[NUMBER_OF_BEAMS*4];
float scat_exp[NUMBER_OF_BEAMS*4];
float scat_tot_toa[NUMBER_OF_BEAMS];
float scat_kpc_ant[NUMBER_OF_BEAMS*4]; 
float scat_kpc_toa[NUMBER_OF_BEAMS*4];
float scat_tot_kpc_toa[NUMBER_OF_BEAMS];
float scat_swind[NUMBER_OF_BEAMS];
float scat_swindunc[NUMBER_OF_BEAMS];
float scat_esurf[NUMBER_OF_BEAMS*2];
float scat_esurf_unc[NUMBER_OF_BEAMS*2];
float scat_sLIfrac[NUMBER_OF_BEAMS*2];
float scat_sCenLon[MAXCYC*NUMBER_OF_BEAMS];
float scat_sCenLat[MAXCYC*NUMBER_OF_BEAMS];
float scat_sEdgeLon[MAXCYC*NUMBER_OF_BEAMS*4];
float scat_sEdgeLat[MAXCYC*NUMBER_OF_BEAMS*4];
float scat_polarization_roll[MAXCYC*NUMBER_OF_BEAMS];

uint32_t scat_sFlags[NUMBER_OF_BEAMS];
uint8_t  scat_sRFI_Flags[NUMBER_OF_BEAMS*2];

uint16_t scat_samples[NUMBER_OF_BEAMS];
uint16_t rad_samples[NUMBER_OF_BEAMS][RADIOMETER_POLARIZATIONS];


// TA_hat is antenna temperatures before loss correction [kelvin] 
// TA     is antenna temperatures after  loss correction [kelvin]
// (1=V, 2=H, 3=3rd Stokes)
float *TA_hat, *TA, *TA0;
TA_hat = new float[MAXCYC*NUMBER_OF_BEAMS*3];
TA     = new float[MAXCYC*NUMBER_OF_BEAMS*3];
TA0    = new float[MAXCYC*NUMBER_OF_BEAMS*3];


// TF_hat is antenna temperatures before loss correction w/i RFI filter  [kelvin] (1=V, 2=H, 3=3rd Stokes)
// TF     is antenna temperatures after  loss correction w/i RFI filter  [kelvin] (1=V, 2=H, 3=3rd Stokes)
float *TF_hat, *TF, *TF0;
TF_hat = new float[MAXCYC*NUMBER_OF_BEAMS*3];
TF     = new float[MAXCYC*NUMBER_OF_BEAMS*3];
TF0    = new float[MAXCYC*NUMBER_OF_BEAMS*3];

float *TA_hat_drift_corr;
TA_hat_drift_corr = new float[MAXCYC*NUMBER_OF_BEAMS*2];

float *cellon, *cellat, *celtht, *celphi, *suntht, *sunphi, *sunglt, *moonglt;
float *glxlon, *glxlat, *cellonfoot, *cellatfoot;
double *zang, *sun_zenith, *sclon, *sclat, *scalt;
double *sund, *sunr, *moond, *moonr, *bore_sight;
uint8_t *acsmode;

cellon = new float[MAXCYC*NUMBER_OF_BEAMS];
cellat = new float[MAXCYC*NUMBER_OF_BEAMS];
celtht = new float[MAXCYC*NUMBER_OF_BEAMS];
celphi = new float[MAXCYC*NUMBER_OF_BEAMS];
suntht = new float[MAXCYC*NUMBER_OF_BEAMS];
sunphi = new float[MAXCYC*NUMBER_OF_BEAMS];
sunglt = new float[MAXCYC*NUMBER_OF_BEAMS];
moonglt = new float[MAXCYC*NUMBER_OF_BEAMS];
glxlon = new float[MAXCYC*NUMBER_OF_BEAMS];
glxlat = new float[MAXCYC*NUMBER_OF_BEAMS];
zang = new double[MAXCYC];
for (iblk=0; iblk<MAXCYC; iblk++) zang[iblk] = -999;
sun_zenith = new double[MAXCYC];
sclon = new double[MAXCYC];
sclat = new double[MAXCYC];
scalt = new double[MAXCYC];

cellonfoot = new float[MAXCYC*NUMBER_OF_BEAMS*4];
cellatfoot = new float[MAXCYC*NUMBER_OF_BEAMS*4];

sund = new double[MAXCYC*NUMBER_OF_BEAMS];
sunr = new double[MAXCYC*NUMBER_OF_BEAMS];
moond = new double[MAXCYC*NUMBER_OF_BEAMS];
moonr = new double[MAXCYC*NUMBER_OF_BEAMS];

bore_sight = new double[MAXCYC*NUMBER_OF_BEAMS*3];

acsmode = new uint8_t[MAXCYC];


float *solar_flux, *solar_xray;
float *fland, *gland, *fice, *gice, *sss_reference, *surtep, *subtep;
float *surp, *cwat, *swh;
float *sm, *swe;
float *winspd_ncep, *winspd_hh, *winspd_hhh, *windir, *tran, *tbup, *tbdw;
float *tb_landcorr_vpol, *tb_landcorr_hpol;

solar_flux = new float[MAXCYC];
solar_xray = new float[MAXCYC];
fland = new float[MAXCYC*NUMBER_OF_BEAMS];
gland = new float[MAXCYC*NUMBER_OF_BEAMS];
fice = new float[MAXCYC*NUMBER_OF_BEAMS];
gice = new float[MAXCYC*NUMBER_OF_BEAMS];
swh = new float[MAXCYC*NUMBER_OF_BEAMS];
sss_reference = new float[MAXCYC*NUMBER_OF_BEAMS];
tb_landcorr_vpol = new float[MAXCYC*NUMBER_OF_BEAMS];
tb_landcorr_hpol = new float[MAXCYC*NUMBER_OF_BEAMS];
surtep = new float[MAXCYC*NUMBER_OF_BEAMS];
subtep = new float[MAXCYC*NUMBER_OF_BEAMS];
surp = new float[MAXCYC*NUMBER_OF_BEAMS];
cwat = new float[MAXCYC*NUMBER_OF_BEAMS];
sm = new float[MAXCYC*NUMBER_OF_BEAMS];
swe = new float[MAXCYC*NUMBER_OF_BEAMS];
winspd_ncep = new float[MAXCYC*NUMBER_OF_BEAMS];
winspd_hh = new float[MAXCYC*NUMBER_OF_BEAMS];
winspd_hhh = new float[MAXCYC*NUMBER_OF_BEAMS];
windir = new float[MAXCYC*NUMBER_OF_BEAMS];
tran = new float[MAXCYC*NUMBER_OF_BEAMS];
tbup = new float[MAXCYC*NUMBER_OF_BEAMS];
tbdw = new float[MAXCYC*NUMBER_OF_BEAMS];


float *tagal_dir, *tagal_ref, *tasun_dir, *tasun_ref, *tasun_bak, *tamon_ref;
tagal_dir = new float[MAXCYC*NUMBER_OF_BEAMS*3];
tagal_ref = new float[MAXCYC*NUMBER_OF_BEAMS*3];
tasun_dir = new float[MAXCYC*NUMBER_OF_BEAMS*3];
tasun_ref = new float[MAXCYC*NUMBER_OF_BEAMS*3];
tasun_bak = new float[MAXCYC*NUMBER_OF_BEAMS*3];
tamon_ref = new float[MAXCYC*NUMBER_OF_BEAMS*3];


///////////////////////////////////////////////


int32_t prodOffset[MAXNAQPROD];

int8_t prodMultiplicity[MAXNAQPROD];

float *prodPtr[MAXNAQPROD];

// Note: If prodMultiplicity is negative then the associated variable 
//       has no "iblk", ie, [MAXCYC dimension.

prodPtr[TAV] = TA;
prodMultiplicity[TAV] = 3;
prodOffset[TAV] = 0;

prodPtr[TAH] = TA;
prodMultiplicity[TAH] = 3;
prodOffset[TAH] = 1;

prodPtr[TA3] = TA;
prodMultiplicity[TA3] = 3;
prodOffset[TA3] = 2;


prodPtr[TFV] = TF;
prodMultiplicity[TFV] = 3;
prodOffset[TFV] = 0;

prodPtr[TFH] = TF;
prodMultiplicity[TFH] = 3;
prodOffset[TFH] = 1;

prodPtr[TF3] = TF;
prodMultiplicity[TF3] = 3;
prodOffset[TF3] = 2;



prodPtr[TAV0] = TA0;
prodMultiplicity[TAV0] = 3;
prodOffset[TAV0] = 0;

prodPtr[TAH0] = TA0;
prodMultiplicity[TAH0] = 3;
prodOffset[TAH0] = 1;

prodPtr[TA30] = TA0;
prodMultiplicity[TA30] = 3;
prodOffset[TA30] = 2;


prodPtr[TFV0] = TF0;
prodMultiplicity[TFV0] = 3;
prodOffset[TFV0] = 0;

prodPtr[TFH0] = TF0;
prodMultiplicity[TFH0] = 3;
prodOffset[TFH0] = 1;

prodPtr[TF30] = TF0;
prodMultiplicity[TF30] = 3;
prodOffset[TF30] = 2;


prodPtr[TOIV] = tb_toi;
prodMultiplicity[TOIV] = -3;
prodOffset[TOIV] = 0;

prodPtr[TOIH] = tb_toi;
prodMultiplicity[TOIH] = -3;
prodOffset[TOIH] = 1;

prodPtr[TOI3] = tb_toi;
prodMultiplicity[TOI3] = -3;
prodOffset[TOI3] = 2;


prodPtr[TOAV] = tb_toa_nolc;
prodMultiplicity[TOAV] = -2;
prodOffset[TOAV] = 0;

prodPtr[TOAH] = tb_toa_nolc;
prodMultiplicity[TOAH] = -2;
prodOffset[TOAH] = 1;


prodPtr[FARTAH] = faraday_deg;
prodMultiplicity[FARTAH] = 1;
prodOffset[FARTAH] = 0;


prodPtr[TAGALDV] = tagal_dir;
prodMultiplicity[TAGALDV] = 3;
prodOffset[TAGALDV] = 0;

prodPtr[TAGALDH] = tagal_dir;
prodMultiplicity[TAGALDH] = 3;
prodOffset[TAGALDH] = 1;

prodPtr[TAGALD3] = tagal_dir;
prodMultiplicity[TAGALD3] = 3;
prodOffset[TAGALD3] = 2;


prodPtr[TAGALRV] = tagal_ref;
prodMultiplicity[TAGALRV] = 3;
prodOffset[TAGALRV] = 0;

prodPtr[TAGALRH] = tagal_ref;
prodMultiplicity[TAGALRH] = 3;
prodOffset[TAGALRH] = 1;

prodPtr[TAGALR3] = tagal_ref;
prodMultiplicity[TAGALR3] = 3;
prodOffset[TAGALR3] = 2;


prodPtr[TAGALRGOV] = tagal_ref_GO;
prodMultiplicity[TAGALRGOV] = -3;
prodOffset[TAGALRGOV] = 0;

prodPtr[TAGALRGOH] = tagal_ref_GO;
prodMultiplicity[TAGALRGOH] = -3;
prodOffset[TAGALRGOH] = 1;


prodPtr[DTAGALV] = dtagal_ref;
prodMultiplicity[DTAGALV] = -3;
prodOffset[DTAGALV] = 0;

prodPtr[DTAGALH] = dtagal_ref;
prodMultiplicity[DTAGALH] = -3;
prodOffset[DTAGALH] = 1;


prodPtr[TASUNDV] = tasun_dir;
prodMultiplicity[TASUNDV] = 3;
prodOffset[TASUNDV] = 0;

prodPtr[TASUNDH] = tasun_dir;
prodMultiplicity[TASUNDH] = 3;
prodOffset[TASUNDH] = 1;

prodPtr[TASUND3] = tasun_dir;
prodMultiplicity[TASUND3] = 3;
prodOffset[TASUND3] = 2;


prodPtr[TASUNRV] = tasun_ref;
prodMultiplicity[TASUNRV] = 3;
prodOffset[TASUNRV] = 0;

prodPtr[TASUNRH] = tasun_ref;
prodMultiplicity[TASUNRH] = 3;
prodOffset[TASUNRH] = 1;

prodPtr[TASUNR3] = tasun_ref;
prodMultiplicity[TASUNR3] = 3;
prodOffset[TASUNR3] = 2;


prodPtr[TASUNBV] = tasun_bak;
prodMultiplicity[TASUNBV] = 3;
prodOffset[TASUNBV] = 0;

prodPtr[TASUNBH] = tasun_bak;
prodMultiplicity[TASUNBH] = 3;
prodOffset[TASUNBH] = 1;

prodPtr[TASUNB3] = tasun_bak;
prodMultiplicity[TASUNB3] = 3;
prodOffset[TASUNB3] = 2;


prodPtr[TAMONRV] = tamon_ref;
prodMultiplicity[TAMONRV] = 3;
prodOffset[TAMONRV] = 0;

prodPtr[TAMONRH] = tamon_ref;
prodMultiplicity[TAMONRH] = 3;
prodOffset[TAMONRH] = 1;

prodPtr[TAMONR3] = tamon_ref;
prodMultiplicity[TAMONR3] = 3;
prodOffset[TAMONR3] = 2;


prodPtr[TBV] = tb_sur_nolc;
prodMultiplicity[TBV] = 2;
prodOffset[TBV] = 0;

prodPtr[TBH] = tb_sur_nolc;
prodMultiplicity[TBH] = 2;
prodOffset[TBH] = 1;

prodPtr[TBVRC] = tb_sur_rc_nolc;
prodMultiplicity[TBVRC] = 2;
prodOffset[TBVRC] = 0;

prodPtr[TBHRC] = tb_sur_rc_nolc;
prodMultiplicity[TBHRC] = 2;
prodOffset[TBHRC] = 1;

prodPtr[TBVRCLC] = tb_sur_rc;
prodMultiplicity[TBVRCLC] = 2;
prodOffset[TBVRCLC] = 0;

prodPtr[TBHRCLC] = tb_sur_rc;
prodMultiplicity[TBHRCLC] = 2;
prodOffset[TBHRCLC] = 1;

prodPtr[TBVNORCLC] = tb_sur;
prodMultiplicity[TBVNORCLC] = 2;
prodOffset[TBVNORCLC] = 0;

prodPtr[TBHNORCLC] = tb_sur;
prodMultiplicity[TBHNORCLC] = 2;
prodOffset[TBHNORCLC] = 1;

prodPtr[TBCON] = tb_con_nolc;
prodMultiplicity[TBCON] = 1;
prodOffset[TBCON] = 0;

prodPtr[TBCONLC] = tb_con;
prodMultiplicity[TBCONLC] = 1;
prodOffset[TBCONLC] = 0;


prodPtr[SSS] = sss_aquarius_nolc;
prodMultiplicity[SSS] = 1;
prodOffset[SSS] = 0;

prodPtr[SSSLC] = sss_aquarius;
prodMultiplicity[SSSLC] = 1;
prodOffset[SSSLC] = 0;

prodPtr[SSSERRRAN] = sss_aquarius_err_ran;
prodMultiplicity[SSSERRRAN] = 1;
prodOffset[SSSERRRAN] = 0;

prodPtr[SSSERRSYS] = sss_aquarius_err_sys;
prodMultiplicity[SSSERRSYS] = 1;
prodOffset[SSSERRSYS] = 0;

prodPtr[SSSADJ] = sss_aquarius_adj;
prodMultiplicity[SSSADJ] = 1;
prodOffset[SSSADJ] = 0;

prodPtr[DENSITY] = density;
prodMultiplicity[DENSITY] = 1;
prodOffset[DENSITY] = 0;

prodPtr[SPICE] = spiciness;
prodMultiplicity[SPICE] = 1;
prodOffset[SPICE] = 0;


prodPtr[VSM] = vol_sm;
prodMultiplicity[VSM] = 1;
prodOffset[VSM] = 0;


prodPtr[TOAVLC] = tb_toa;
prodMultiplicity[TOAVLC] = -2;
prodOffset[TOAVLC] = 0;

prodPtr[TOAHLC] = tb_toa;
prodMultiplicity[TOAHLC] = -2;
prodOffset[TOAHLC] = 1;

//prodPtr[DSSSLC] = d_sss_aquarius_land_corr;
//prodMultiplicity[DSSSLC] = 1;
//prodOffset[DSSSLC] = 0;


prodPtr[EXPTAV] = ta_expected;
prodMultiplicity[EXPTAV] = 3;
prodOffset[EXPTAV] = 0;

prodPtr[EXPTAH] = ta_expected;
prodMultiplicity[EXPTAH] = 3;
prodOffset[EXPTAH] = 1;

prodPtr[EXPTA3] = ta_expected;
prodMultiplicity[EXPTA3] = 3;
prodOffset[EXPTA3] = 2;


prodPtr[EXPTAV_HHH] = ta_expected_hhh;
prodMultiplicity[EXPTAV_HHH] = 3;
prodOffset[EXPTAV_HHH] = 0;

prodPtr[EXPTAH_HHH] = ta_expected_hhh;
prodMultiplicity[EXPTAH_HHH] = 3;
prodOffset[EXPTAH_HHH] = 1;

prodPtr[EXPTA3_HHH] = ta_expected_hhh;
prodMultiplicity[EXPTA3_HHH] = 3;
prodOffset[EXPTA3_HHH] = 2;


prodPtr[EXPTBV] = tb_expected;
prodMultiplicity[EXPTBV] = 2;
prodOffset[EXPTBV] = 0;

prodPtr[EXPTBH] = tb_expected;
prodMultiplicity[EXPTBH] = 2;
prodOffset[EXPTBH] = 1;

prodPtr[EXPTBV0] = tb_expected0;
prodMultiplicity[EXPTBV0] = 2;
prodOffset[EXPTBV0] = 0;

prodPtr[EXPTBH0] = tb_expected0;
prodMultiplicity[EXPTBH0] = 2;
prodOffset[EXPTBH0] = 1;


prodPtr[HHWINDSPD] = winspd_hh;
prodMultiplicity[HHWINDSPD] = 1;
prodOffset[HHWINDSPD] = 0;

prodPtr[HHHWINDSPD] = winspd_hhh;
prodMultiplicity[HHHWINDSPD] = 1;
prodOffset[HHHWINDSPD] = 0;


prodPtr[DTBSSTWSPDV] = delta_tb_sst_wspd;
prodMultiplicity[DTBSSTWSPDV] = 2;
prodOffset[DTBSSTWSPDV] = 0;

prodPtr[DTBSSTWSPDH] = delta_tb_sst_wspd;
prodMultiplicity[DTBSSTWSPDH] = 2;
prodOffset[DTBSSTWSPDH] = 1;


prodPtr[HVANT] = scat_ant;
prodMultiplicity[HVANT] = -4;
prodOffset[HVANT] = 0;

prodPtr[VVANT] = scat_ant;
prodMultiplicity[VVANT] = -4;
prodOffset[VVANT] = 1;

prodPtr[VHANT] = scat_ant;
prodMultiplicity[VHANT] = -4;
prodOffset[VHANT] = 2;

prodPtr[HHANT] = scat_ant;
prodMultiplicity[HHANT] = -4;
prodOffset[HHANT] = 3;



prodPtr[HVTOA] = scat_toa;
prodMultiplicity[HVTOA] = -4;
prodOffset[HVTOA] = 0;

prodPtr[VVTOA] = scat_toa;
prodMultiplicity[VVTOA] = -4;
prodOffset[VVTOA] = 1;

prodPtr[VHTOA] = scat_toa;
prodMultiplicity[VHTOA] = -4;
prodOffset[VHTOA] = 2;

prodPtr[HHTOA] = scat_toa;
prodMultiplicity[HHTOA] = -4;
prodOffset[HHTOA] = 3;


prodPtr[HVEXP] = scat_exp;
prodMultiplicity[HVEXP] = -4;
prodOffset[HVEXP] = 0;

prodPtr[VVEXP] = scat_exp;
prodMultiplicity[VVEXP] = -4;
prodOffset[VVEXP] = 1;

prodPtr[VHEXP] = scat_exp;
prodMultiplicity[VHEXP] = -4;
prodOffset[VHEXP] = 2;

prodPtr[HHEXP] = scat_exp;
prodMultiplicity[HHEXP] = -4;
prodOffset[HHEXP] = 3;


prodPtr[TOTTOA] = scat_tot_toa;
prodMultiplicity[TOTTOA] = -1;
prodOffset[TOTTOA] = 0;


prodPtr[HVKPCANT] = scat_kpc_ant;
prodMultiplicity[HVKPCANT] = -4;
prodOffset[HVKPCANT] = 0;

prodPtr[VVKPCANT] = scat_kpc_ant;
prodMultiplicity[VVKPCANT] = -4;
prodOffset[VVKPCANT] = 1;

prodPtr[VHKPCANT] = scat_kpc_ant;
prodMultiplicity[VHKPCANT] = -4;
prodOffset[VHKPCANT] = 2;

prodPtr[HHKPCANT] = scat_kpc_ant;
prodMultiplicity[HHKPCANT] = -4;
prodOffset[HHKPCANT] = 3;



prodPtr[HVKPCTOA] = scat_kpc_toa;
prodMultiplicity[HVKPCTOA] = -4;
prodOffset[HVKPCTOA] = 0;

prodPtr[VVKPCTOA] = scat_kpc_toa;
prodMultiplicity[VVKPCTOA] = -4;
prodOffset[VVKPCTOA] = 1;

prodPtr[VHKPCTOA] = scat_kpc_toa;
prodMultiplicity[VHKPCTOA] = -4;
prodOffset[VHKPCTOA] = 2;

prodPtr[HHKPCTOA] = scat_kpc_toa;
prodMultiplicity[HHKPCTOA] = -4;
prodOffset[HHKPCTOA] = 3;


prodPtr[TOTKPCTOA] = scat_tot_kpc_toa;
prodMultiplicity[TOTKPCTOA] = -1;
prodOffset[TOTKPCTOA] = 0;


prodPtr[SWINDSPD] = scat_swind;
prodMultiplicity[SWINDSPD] = -1;
prodOffset[SWINDSPD] = 0;

prodPtr[SWINDUNC] = scat_swindunc;
prodMultiplicity[SWINDUNC] = -1;
prodOffset[SWINDUNC] = 0;


prodPtr[ESURFV] = scat_esurf;
prodMultiplicity[ESURFV] = -2;
prodOffset[ESURFV] = 0;

prodPtr[ESURFH] = scat_esurf;
prodMultiplicity[ESURFH] = -2;
prodOffset[ESURFH] = 1;

prodPtr[ESURFVUNC] = scat_esurf_unc;
prodMultiplicity[ESURFVUNC] = -2;
prodOffset[ESURFVUNC] = 0;

prodPtr[ESURFHUNC] = scat_esurf_unc;
prodMultiplicity[ESURFHUNC] = -2;
prodOffset[ESURFHUNC] = 1;


prodPtr[AWINDSPD] = winspd_ncep;
prodMultiplicity[AWINDSPD] =  1;
prodOffset[AWINDSPD] = 0;

prodPtr[AWINDDIR] = windir;
prodMultiplicity[AWINDDIR] = 1;
prodOffset[AWINDDIR] = 0;


prodPtr[ASURTEMP] = surtep;
prodMultiplicity[ASURTEMP] = 1;
prodOffset[ASURTEMP] = 0;

prodPtr[ASURP] = surp;
prodMultiplicity[ASURP] = 1;
prodOffset[ASURP] = 0;

prodPtr[ACWAT] = cwat;
prodMultiplicity[ACWAT] = 1;
prodOffset[ACWAT] = 0;

prodPtr[ASM] = sm;
prodMultiplicity[ASM] = 1;
prodOffset[ASM] = 0;

prodPtr[ASUBTEMP] = subtep;
prodMultiplicity[ASUBTEMP] = 1;
prodOffset[ASUBTEMP] = 0;

prodPtr[ASWE] = swe;
prodMultiplicity[ASWE] = 1;
prodOffset[ASWE] = 0;

prodPtr[ASSS] = sss_reference;
prodMultiplicity[ASSS] = 1;
prodOffset[ASSS] = 0;


prodPtr[ATRANS] = tran;
prodMultiplicity[ATRANS] = 1;
prodOffset[ATRANS] = 0;

prodPtr[ATBUP] = tbup;
prodMultiplicity[ATBUP] = 1;
prodOffset[ATBUP] = 0;

prodPtr[ATBDW] = tbdw;
prodMultiplicity[ATBDW] = 1;
prodOffset[ATBDW] = 0;

prodPtr[ASWH] = swh;
prodMultiplicity[ASWH] = 1;
prodOffset[ASWH] = 0;


prodPtr[RLANDFRC] = gland;
prodMultiplicity[RLANDFRC] = 1;
prodOffset[RLANDFRC] = 0;

prodPtr[RICEFRC] = gice;
prodMultiplicity[RICEFRC] = 1;
prodOffset[RICEFRC] = 0;

//prodPtr[RGICE] = gice;
//prodMultiplicity[RGICE] = 1;
//prodOffset[RGICE] = 0;


prodPtr[SLANDFRC] = scat_sLIfrac;
prodMultiplicity[SLANDFRC] = -2;
prodOffset[RLANDFRC] = 0;

prodPtr[SICEFRC] = scat_sLIfrac;
prodMultiplicity[SICEFRC] = -2;
prodOffset[SICEFRC] = 1;

