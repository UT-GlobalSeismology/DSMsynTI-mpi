/*************************************************************************
 *  ************** spcsac ****************
 *  Tool to transform frequency spectrums computed by the DSM software
 *  into time domain data in a SAC binary format.
 *
 *  Main historical authors: N.Fuji, K.Kawai, N.Takeuchi, R.J.Geller
 *  (C) 2006.3  University of Tokyo
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation, either version 3 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program.  If not, see <https://www.gnu.org/licenses/>.
 *
 **************************************************************************/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>
#define PI 3.14159265358979

struct sacheader {
    double DELTA;
    double DEPMIN;
    double DEPMAX;
    double SCALE;
    double ODELTA;
    double B;
    double E;
    double O;
    double A;
    double H9;
    double T0;
    double T1;
    double T2;
    double T3;
    double T4;
    double T5;
    double T6;
    double T7;
    double T8;
    double T9;
    double F;
    double RESP0;
    double RESP1;
    double RESP2;
    double RESP3;
    double RESP4;
    double RESP5;
    double RESP6;
    double RESP7;
    double RESP8;
    double RESP9;
    double STLA;
    double STLO;
    double STEL;
    double STDP;
    double EVLA;
    double EVLO;
    double EVEL;
    double EVDP;
    double H39;
    double USER0;
    double USER1;
    double USER2;
    double USER3;
    double USER4;
    double USER5;
    double USER6;
    double USER7;
    double USER8;
    double USER9;
    double DIST;
    double AZ;
    double BAZ;
    double GCARC;
    double H54;
    double H55;
    double DEPMEN;
    double CMPAZ;
    double CMPINC;
    double H59;
    double H60;
    double H61;
    double H62;
    double H63;
    double H64;
    double H65;
    double H66;
    double H67;
    double H68;
    double H69;
    long NZYEAR;
    long NZJDAY;
    long NZHOUR;
    long NZMIN;
    long NZSEC;
    long NZMSEC;
    long NVHDR;
    long H77;
    long H78;
    long NPTS;
    long H80;
    long H81;
    long H82;
    long H83;
    long H84;
    long IFTYPE;
    long IDEP;
    long IZTYPE;
    long H88;
    long H89;
    long ISTREG;
    long IEVREG;
    long IEVTYP;
    long IQUAL;
    long ISYNTH;
    long H95;
    long H96;
    long H97;
    long H98;
    long H99;
    long H100;
    long H101;
    long H102;
    long H103;
    long H104;
    long LEVEN;
    long LPSPOL;
    long LOVROK;
    long LCALDA;
    long H109;
    char KSTNM[9];
    char KEVNM[17];
    char KHOLE[9];
    char KO[9];
    char KA[9];
    char KT0[9];
    char KT1[9];
    char KT2[9];
    char KT3[9];
    char KT4[9];
    char KT5[9];
    char KT6[9];
    char KT7[9];
    char KT8[9];
    char KT9[9];
    char KF[9];
    char KUSER0[9];
    char KUSER1[9];
    char KUSER2[9];
    char KCMPNM[9];
    char KNETWK[9];
    char KDATRD[9];
    char KINST[9];
};

void onespcsac(float samplingFreq, int icomplexinv, double muRperturb, double Qperturb, int psvorsh, int component, long lsmooth,
               char *backr, char *backt, char *backz, char *backp, char *backs,
               char *psvfile, char *shfile, char *spcfile, int asciioutnput);
void gcc2ggc(double *theta);
void helpmessage(void);
void lsmoothfinder(void);
void nojokenoscience(void);
int filexist(const char *filename);
void error01(char *spcfile);
void modeprompt(int *filordir, int *psvorsh, int *component,
                char *backr, char *backt, char *backz, char *backp, char *backs, char *spcfile);
void cdft(long n, int isgn, double *a, int *ip, double *w);
void rdft(long n, int isgn, double *a, int *ip, double *w);
void ddct(long n, int isgn, double *a, int *ip, double *w);
void ddst(long n, int isgn, double *a, int *ip, double *w);
void dfct(long n, double *a, double *t, int *ip, double *w);
void dfst(long n, double *a, double *t, int *ip, double *w);
void makewt(long nw, int *ip, double *w);
void makect(long nc, int *ip, double *c);
void bitrv2(long n, int *ip, double *a);
void bitrv2conj(long n, int *ip, double *a);
void cftfsub(long n, double *a, double *w);
void cftbsub(long n, double *a, double *w);
void cft1st(long n, double *a, double *w);
void cftmdl(long n, int l, double *a, double *w);
void rftfsub(long n, double *a, int nc, double *c);
void rftbsub(long n, double *a, int nc, double *c);
void dctsub(long n, double *a, int nc, double *c);
void dstsub(long n, double *a, int nc, double *c);
int lsmoothfinder_(float tlen, int np0, float freq);
double *readsac(char *sacfile, struct sacheader *hv, long nptsmodifier);
int writesac(char *sacfile, struct sacheader *hv, double *sacdata);
int writeasciisac(char *sacfile, struct sacheader *hv, double *sacdata);
int newsacheader(struct sacheader *hv);

int main(int argc, char *argv[]) {
    int ch;
    extern char *optarg;
    extern int optind, oterr;
    FILE *file_spc;
    char spcfile[80], psvfile[80], shfile[80];
    char command[100];
    int filordir = 0;
    int psvorsh = 3;
    int component = 7;
    long lsmooth = 4;

    int icomplexinv = 0;
    double muRperturb, Qperturb;
    int i;
    float samplingFreq = 0.0;

    char backr[20];
    char backt[20];
    char backz[20];
    char backp[20];
    char backs[20];

    int asciioutput = 0;

    strcpy(psvfile, "\0");
    strcpy(shfile, "\0");

    strcpy(backr, ".Rs");
    strcpy(backt, ".Ts");
    strcpy(backz, ".Zs");
    strcpy(backp, ".PSV.spc");
    strcpy(backs, ".SH.spc");
    strcpy(spcfile, "./");

    while ((ch = getopt(argc, argv, "M:Q:H:hevid:f:l:m:c:r:t:z:p:s:a")) != -1) {
        switch (ch) {
            case 'Q':
                icomplexinv = 1;
                Qperturb = atof(optarg);
                printf("Q=%lf\n", Qperturb);
                break;
            case 'M':
                icomplexinv = 1;
                muRperturb = atof(optarg);
                printf("mu=%lf\n", muRperturb);
                break;
            case 'H':
                samplingFreq = atof(optarg);
                printf("sampling Hz = %lf\n", samplingFreq);
                break;
            case 'h':
                helpmessage();
                break;
            case 'e':
                modeprompt(&filordir, &psvorsh, &component, backr, backt, backz, backp, backs, spcfile);
                break;
            case 'v':
                puts("<< spcsac, Release 1.0.0 >>");
                nojokenoscience();
                break;
            case 'i':
                lsmoothfinder();
                break;
            case 'a':
                asciioutput = 1;
                break;
            case 'd':
                strcpy(spcfile, optarg);
                break;
            case 'f':
                filordir = 1;
                strcpy(spcfile, optarg);
                strcpy(backr, "tmpsac.Rs");
                strcpy(backt, "tmpsac.Ts");
                strcpy(backz, "tmpsac.Zs");
                break;
            case 'l':
                lsmooth = atol(optarg);
                break;
            case 'm':
                if (!strcmp(optarg, "PSV")) psvorsh = 1;
                if (!strcmp(optarg, "SH")) psvorsh = 2;
                if (!strcmp(optarg, "psv")) psvorsh = 1;
                if (!strcmp(optarg, "sh")) psvorsh = 2;
                break;
            case 'c':
                if (!strcmp(optarg, "T")) component = 4;
                if (!strcmp(optarg, "R")) component = 2;
                if (!strcmp(optarg, "Z")) component = 1;
                if (!strcmp(optarg, "TR")) component = 6;
                if (!strcmp(optarg, "RT")) component = 6;
                if (!strcmp(optarg, "TZ")) component = 5;
                if (!strcmp(optarg, "ZT")) component = 5;
                if (!strcmp(optarg, "RZ")) component = 3;
                if (!strcmp(optarg, "ZR")) component = 3;
                break;
            case 'r':
                strcpy(backr, optarg);
                break;
            case 't':
                strcpy(backt, optarg);
                break;
            case 'z':
                strcpy(backz, optarg);
                break;
            case 'p':
                strcpy(backp, optarg);
                break;
            case 's':
                strcpy(backs, optarg);
                break;
            default:
                break;
        }
        optarg = NULL;
    }

    argc -= optind;
    argv += optind;

    if (argc != 0) helpmessage();

    if (filordir == 1) {
        psvorsh = 1;
        strcpy(psvfile, spcfile);
        strcpy(spcfile, "\0");
        onespcsac(samplingFreq, icomplexinv, muRperturb, Qperturb, psvorsh, component, lsmooth, backr, backt, backz, backp, backs,
                  psvfile, shfile, spcfile, asciioutput);
    }

    if (filordir == 0) {
        sprintf(command, "cd %s", spcfile);
        if (system(command)) {
            printf("Somehow we cannot get the directory %s.\n", spcfile);
            exit(1);
        }
        if ((psvorsh == 3) || (psvorsh == 1)) {
            if (filexist("file_spc")) system("rm file_spc");
            sprintf(command, "for file in %s/*%s; \n do \n echo $file >>file_spc; \n done\n", spcfile, backp);
            if (system(command)) {
                printf("Somehow we cannot find spc files %s/*%s.\n", spcfile, backp);
                exit(1);
            }
            if (!(file_spc = fopen("file_spc", "r"))) {
                puts("Somehow file_spc was not produced.");
                exit(1);
            }
        }
        if (psvorsh == 2) {
            if (filexist("file_spc")) system("rm file_spc");
            sprintf(command, "for file in %s/*%s; \n do \n echo $file >>file_spc; \n done\n", spcfile, backs);
            if (system(command)) {
                printf("Somehow we cannot find spc files %s/*%s.\n", spcfile, backs);
                exit(1);
            }
            if (!(file_spc = fopen("file_spc", "r"))) {
                puts("Somehow file_spc was not produced.");
                exit(1);
            }
        }
        strcpy(spcfile, "\0");
        while (fgets(spcfile, 80, file_spc) != NULL) {
            if (psvorsh == 3) {
                spcfile[strlen(spcfile) - strlen(backp) - 1] = '\0';
                sprintf(psvfile, "%s%s", spcfile, backp);
                sprintf(shfile, "%s%s", spcfile, backs);
            }
            if (psvorsh == 2) {
                spcfile[strlen(spcfile) - strlen(backs) - 1] = '\0';
                sprintf(shfile, "%s%s", spcfile, backs);
            }
            if (psvorsh == 1) {
                spcfile[strlen(spcfile) - strlen(backp) - 1] = '\0';
                sprintf(psvfile, "%s%s", spcfile, backp);
            }
            onespcsac(samplingFreq, icomplexinv, muRperturb, Qperturb, psvorsh, component, lsmooth, backr, backt, backz, backp, backs,
                      psvfile, shfile, spcfile, asciioutput);
        }
    }

    if (filexist("file_spc")) system("rm file_spc");

    puts(" ----------------------------------------------------------------------");
    puts("  ARIGATO, VIELEN DANK, XIEXIE, MERCI, OBRIGADO, DOBRE, THANK YOU...   ");

    return 0;
}

void onespcsac(float samplingFreq, int icomplexinv, double muRperturb, double Qperturb, int psvorsh, int component, long lsmooth,
               char *backr, char *backt, char *backz, char *backp, char *backs,
               char *psvfile, char *shfile, char *spcfile, int asciioutput) {
    FILE *file_psv, *file_sh;
    double tlen;
    long np0, np = 1;
    long i;
    double omegai, theta, phi;
    long nbody, ncomp;

    double *w;
    int *ip;
    double *x;
    double ftmp;
    long nerr, n1, m1;

    double *psvr, *psvt, *psvp;
    double *shr, *sht, *shp;
    double *ypsvr, *yshr, *ynetr;
    double *ypsvt, *ysht, *ynett;
    double *ypsvp, *yshp, *ynetp;
    char output[80];
    double stla, stlo, evla, evlo, r0;
    long tmpint;
    double tmpDouble0, tmpDouble1;
    struct sacheader hv;

    evla = 89.99999;
    evlo = 0.0;

    r0 = 123456789;

    if ((psvorsh == 3) || (psvorsh == 1)) {
        file_psv = fopen(psvfile, "r");
        if (!fscanf(file_psv, "%lf\n", &tlen)) error01(psvfile);
        fscanf(file_psv, "%ld %ld %ld\n", &np0, &nbody, &ncomp);
        fscanf(file_psv, "%lf %lf %lf\n", &omegai, &stla, &stlo);
        fscanf(file_psv, "%lf %lf %lf\n", &evla, &evlo, &r0);
    }
    if ((psvorsh == 3) || (psvorsh == 2)) {
        file_sh = fopen(shfile, "r");
        if (!fscanf(file_sh, "%lf\n", &tlen)) error01(shfile);
        fscanf(file_sh, "%ld %ld %ld\n", &np0, &nbody, &ncomp);
        fscanf(file_sh, "%lf %lf %lf\n", &omegai, &stla, &stlo);
        fscanf(file_sh, "%lf %lf %lf\n", &evla, &evlo, &r0);
    }

    r0 = 6371.0 - r0;

    if (samplingFreq != 0.0) {
        lsmooth = lsmoothfinder_(tlen, np0, samplingFreq);
    }

    i = 1;
    while (i < lsmooth) i = i * 2;
    lsmooth = i;
    i = 0;

    np = 1;
    while (np < np0) np = np * 2;
    np = np * lsmooth;

    if (psvorsh == 3) {
        printf("spcsac is about to transform %s and %s.\n", psvfile, shfile);
        printf("sampling frequency is %f Hz!\n", (2 * np / tlen));
    }

    if (psvorsh == 2) {
        printf("spcsac is about to transform %s.\n", shfile);
        printf("sampling frequency is %f Hz!\n", (2 * np / tlen));
    }

    if (psvorsh == 1) {
        printf("(spcsac is about to transform %s.\n", psvfile);
        printf("sampling frequency is %f Hz!\n", (2 * np / tlen));
    }

    if ((psvorsh == 3) || (psvorsh == 1)) {
        psvr = (double *)calloc((size_t)4 * np, sizeof(double));
        if (psvr == NULL) {
            printf("We could not obtain enough memory for psvr.\n");
            exit(1);
        }
        psvt = (double *)calloc((size_t)4 * np, sizeof(double));
        if (psvt == NULL) {
            printf("We could not obtain enough memory for psvt.\n");
            exit(1);
        }
        psvp = (double *)calloc((size_t)4 * np, sizeof(double));
        if (psvp == NULL) {
            printf("We could not obtain enough memory for psvp.\n");
            exit(1);
        }
        if (icomplexinv == 0) {
            for (i = 0; i <= np0; i++) {
                if (!fscanf(file_psv, "%ld %lf %lf\n", &tmpint, &tmpDouble0, &tmpDouble1))
                    error01(psvfile);
                psvr[2 * tmpint] = tmpDouble0;
                psvr[2 * tmpint + 1] = tmpDouble1;
                if (!fscanf(file_psv, "%lf %lf\n", &tmpDouble0, &tmpDouble1))
                    error01(psvfile);
                psvt[2 * tmpint] = tmpDouble0;
                psvt[2 * tmpint + 1] = tmpDouble1;
                if (!fscanf(file_psv, "%lf %lf\n", &tmpDouble0, &tmpDouble1))
                    error01(psvfile);
                psvp[2 * tmpint] = tmpDouble0;
                psvp[2 * tmpint + 1] = tmpDouble1;
            }
        } else if (icomplexinv == 1) {
            for (i = 0; i <= np0; i++) {
                if (!fscanf(file_psv, "%ld %lf %lf\n", &tmpint, &tmpDouble0, &tmpDouble1))
                    error01(psvfile);
                psvr[2 * tmpint] = tmpDouble0;
                psvr[2 * tmpint + 1] = tmpDouble1;
                if (!fscanf(file_psv, "%lf %lf\n", &tmpDouble0, &tmpDouble1))
                    error01(psvfile);
                psvt[2 * tmpint] = tmpDouble0;
                psvt[2 * tmpint + 1] = tmpDouble1;
                if (!fscanf(file_psv, "%lf %lf\n", &tmpDouble0, &tmpDouble1))
                    error01(psvfile);
                psvp[2 * tmpint] = tmpDouble0;
                psvp[2 * tmpint + 1] = tmpDouble1;

                if (tmpint != 0) {
                    double parreal, parimag, frequency;
                    double parQr, parQi;
                    double factor;
                    frequency = (double)tmpint / tlen;
                    parQr = muRperturb * 2 * log(frequency) / PI;
                    parQi = muRperturb * (1 + 4 * log(frequency) / PI / Qperturb);
                    factor = 1 + 2 * log(frequency) / PI / Qperturb;
                    factor *= sqrt(1 + 1 / (Qperturb * Qperturb));
                    printf("%lf\n", factor);
                    parreal = psvr[2 * tmpint];
                    parimag = psvr[2 * tmpint + 1];
                    psvr[2 * tmpint] = parreal * parQr - parimag * parQi;
                    psvr[2 * tmpint + 1] = parreal * parQi + parimag * parQr;
                    psvr[2 * tmpint] *= 0.001 / factor;
                    psvr[2 * tmpint + 1] *= 0.001 / factor;

                    parreal = psvt[2 * tmpint];
                    parimag = psvr[2 * tmpint + 1];
                    psvt[2 * tmpint] = parreal * parQr - parimag * parQi;
                    psvt[2 * tmpint + 1] = parreal * parQi + parimag * parQr;
                    psvt[2 * tmpint] *= 0.001 / factor;
                    psvt[2 * tmpint] *= 0.001 / factor;

                    parreal = psvp[2 * tmpint];
                    parimag = psvp[2 * tmpint + 1];
                    psvp[2 * tmpint] = parreal * parQr - parimag * parQi;
                    psvp[2 * tmpint + 1] = parreal * parQi + parimag * parQr;
                    psvp[2 * tmpint] *= 0.001 / factor;
                    psvp[2 * tmpint + 1] *= 0.001 / factor;
                }
            }
        }
        fclose(file_psv);
        for (i = 1; i < np; i++) {
            n1 = np + i;
            m1 = np - i;
            psvr[2 * n1] = psvr[2 * m1];
            psvr[2 * n1 + 1] = -psvr[2 * m1 + 1];
            psvt[2 * n1] = psvt[2 * m1];
            psvt[2 * n1 + 1] = -psvt[2 * m1 + 1];
            psvp[2 * n1] = psvp[2 * m1];
            psvp[2 * n1 + 1] = -psvp[2 * m1 + 1];
        }
    }
    if ((psvorsh == 3) || (psvorsh == 2)) {
        shr = (double *)calloc((size_t)4 * np, sizeof(double));
        if (shr == NULL) {
            printf("We could not obtain enough memory for shr.\n");
            exit(1);
        }
        sht = (double *)calloc((size_t)4 * np, sizeof(double));
        if (sht == NULL) {
            printf("We could not obtain enough memory for sht.\n");
            exit(1);
        }
        shp = (double *)calloc((size_t)4 * np, sizeof(double));
        if (shp == NULL) {
            printf("We could not obtain enough memory for shp.\n");
            exit(1);
        }

        if (icomplexinv == 0) {
            for (i = 0; i <= np0; i++) {
                if (!fscanf(file_sh, "%ld %lf %lf\n", &tmpint, &tmpDouble0, &tmpDouble1))
                    error01(shfile);
                shr[2 * tmpint] = tmpDouble0;
                shr[2 * tmpint + 1] = tmpDouble1;
                if (!fscanf(file_sh, "%lf %lf\n", &tmpDouble0, &tmpDouble1))
                    error01(shfile);
                sht[2 * tmpint] = tmpDouble0;
                sht[2 * tmpint + 1] = tmpDouble1;
                if (!fscanf(file_sh, "%lf %lf\n", &tmpDouble0, &tmpDouble1))
                    error01(shfile);
                shp[2 * tmpint] = tmpDouble0;
                shp[2 * tmpint + 1] = tmpDouble1;
            }
        } else if (icomplexinv == 1) {
            /*printf ("%lf, %lf\n", muRperturb, Qperturb);*/

            for (i = 0; i <= np0; i++) {
                if (!fscanf(file_sh, "%ld %lf %lf\n", &tmpint, &tmpDouble0, &tmpDouble1))
                    error01(shfile);
                shr[2 * tmpint] = tmpDouble0;
                shr[2 * tmpint + 1] = tmpDouble1;
                if (!fscanf(file_sh, "%lf %lf\n", &tmpDouble0, &tmpDouble1))
                    error01(shfile);
                sht[2 * tmpint] = tmpDouble0;
                sht[2 * tmpint + 1] = tmpDouble1;
                if (!fscanf(file_sh, "%lf %lf\n", &tmpDouble0, &tmpDouble1))
                    error01(shfile);
                shp[2 * tmpint] = tmpDouble0;
                shp[2 * tmpint + 1] = tmpDouble1;

                if (tmpint != 0) {
                    double parreal, parimag, frequency;
                    double parQr, parQi;
                    double factor;
                    frequency = (double)tmpint / tlen;
                    parQr = muRperturb * 2 * log(frequency) / PI;
                    parQi = muRperturb * (1 + 4 * log(frequency) / PI / Qperturb);

                    factor = 1 + 2 * log(frequency) / PI / Qperturb;
                    factor *= sqrt(1 + 1 / (Qperturb * Qperturb));
                    /*printf("%e %e %e %lf\n", parQr, parQi, Qperturb, factor);
                     */
                    parreal = shr[2 * tmpint];
                    parimag = shr[2 * tmpint + 1];
                    shr[2 * tmpint] = parreal * parQr - parimag * parQi;
                    shr[2 * tmpint + 1] = parreal * parQi + parimag * parQr;
                    shr[2 * tmpint] *= 0.001 / factor;
                    shr[2 * tmpint + 1] *= 0.001 / factor;

                    parreal = sht[2 * tmpint];
                    parimag = sht[2 * tmpint + 1];

                    sht[2 * tmpint] = parreal * parQr - parimag * parQi;
                    sht[2 * tmpint + 1] = parreal * parQi + parimag * parQr;
                    sht[2 * tmpint] *= 0.001 / factor;
                    sht[2 * tmpint + 1] *= 0.001 / factor;

                    parreal = shp[2 * tmpint];
                    parimag = shp[2 * tmpint + 1];
                    shp[2 * tmpint] = parreal * parQr - parimag * parQi;
                    shp[2 * tmpint + 1] = parreal * parQi + parimag * parQr;
                    shp[2 * tmpint] *= 0.001 / factor;
                    shp[2 * tmpint + 1] *= 0.001 / factor;
                }
            }
        }

        fclose(file_sh);

        for (i = 1; i < np; i++) {
            n1 = np + i;
            m1 = np - i;
            shr[2 * n1] = shr[2 * m1];
            shr[2 * n1 + 1] = -shr[2 * m1 + 1];
            sht[2 * n1] = sht[2 * m1];
            sht[2 * n1 + 1] = -sht[2 * m1 + 1];
            shp[2 * n1] = shp[2 * m1];
            shp[2 * n1 + 1] = -shp[2 * m1 + 1];
        }
    }

    np = np * 2;

    if ((component == 7) || (component == 5) || (component == 3) || (component == 1)) {
        w = (double *)calloc((size_t)np * 5 / 2, sizeof(double));
        if (w == NULL) {
            printf("We could not obtain enough memory for w.\n");
            exit(1);
        }

        ip = (int *)calloc((size_t)(2 * sqrt(np) + 4), sizeof(int));
        if (ip == NULL) {
            printf("We could not obtain enough memory for ip.\n");
            exit(1);
        }

        if ((psvorsh == 3) || (psvorsh == 1)) {
            ip[0] = 0;
            cdft(2 * np, 1, psvr, ip, w);
        }
        if ((psvorsh == 3) || (psvorsh == 2)) {
            ip[0] = 0;
            cdft(2 * np, 1, shr, ip, w);
        }

        free(ip);
        free(w);

        x = (double *)calloc((size_t)np, sizeof(double));
        if (x == NULL) {
            printf("We could not obtain enough memory for x.\n");
            exit(1);
        }

        ypsvr = (double *)calloc((size_t)np, sizeof(double));
        if (ypsvr == NULL) {
            printf("We could not obtain enough memory for ypsvr.\n");
            exit(1);
        }

        yshr = (double *)calloc((size_t)np, sizeof(double));
        if (yshr == NULL) {
            printf("We could not obtain enough memory for yshr.\n");
            exit(1);
        }

        ynetr = (double *)calloc((size_t)np, sizeof(double));
        if (ynetr == NULL) {
            printf("We could not obtain enough memory for ynetr.\n");
            exit(1);
        }

        for (i = 0; i < np; i++) {
            x[i] = (double)(tlen * (double)i / (double)np);
            if (psvorsh != 2)
                ypsvr[i] = (double)(psvr[2 * i] * exp(omegai * (double)x[i]) / tlen * 1000.0);
            if (psvorsh != 1)
                yshr[i] = (double)(shr[2 * i] * exp(omegai * (double)x[i]) / tlen * 1000.0);
            ynetr[i] = ypsvr[i] + yshr[i];
        }
        if (psvorsh != 2) free(psvr);
        if (psvorsh != 1) free(shr);

        newsacheader(&hv);
        hv.NPTS = np;
        hv.B = x[0];
        hv.E = x[np - 1];
        hv.DELTA = (double)tlen / np;
        hv.STLA = stla;
        hv.STLO = stlo;
        hv.EVLA = evla;
        hv.EVLO = evlo;
        hv.EVDP = r0;
        strcpy(hv.KCMPNM, "vertical");

        /*
        newhdr();
        setnhv("NPTS",&np,&nerr,4);
        setfhv("B",&x[0],&nerr,1);
        setfhv("E",&x[np-1],&nerr,1);
        ftmp=(float)tlen/np;
        setfhv("DELTA",&ftmp,&nerr,5);
        ftmp=stla;
        setfhv("STLA",&ftmp,&nerr,4);
        ftmp=stlo;
        setfhv("STLO",&ftmp,&nerr,4);
        ftmp=evla;
        setfhv("EVLA",&ftmp,&nerr,4);
        ftmp=evlo;
        setfhv("EVLO",&ftmp,&nerr,4);
        ftmp=r0;
        setfhv("EVDP",&ftmp,&nerr,4);
        */

        if (psvorsh == 3) {
            sprintf(output, "%s%s", spcfile, backz);
            if (asciioutput == 0) {
                writesac(output, &hv, ynetr);
            } else if (asciioutput == 1) {
                writeasciisac(output, &hv, ynetr);
            }
        }
        if (psvorsh == 1) {
            sprintf(output, "%s%s", spcfile, backz);
            if (asciioutput == 0) {
                writesac(output, &hv, ypsvr);
            } else if (asciioutput == 1) {
                writeasciisac(output, &hv, ypsvr);
            }
        }
        if (psvorsh == 2) {
            sprintf(output, "%s%s", spcfile, backz);
            if (asciioutput == 0) {
                writesac(output, &hv, yshr);
            } else if (asciioutput == 1) {
                writeasciisac(output, &hv, yshr);
            }
        }

        free(x);
        free(ypsvr);
        free(yshr);
        free(ynetr);
    } else {
        if (psvorsh != 2) free(psvr);
        if (psvorsh != 1) free(shr);
    }

    if ((component == 7) || (component == 6) || (component == 3) || (component == 2)) {
        w = (double *)calloc((size_t)np * 5 / 2, sizeof(double));
        if (w == NULL) {
            printf("We could not obtain enough memory for w.\n");
            exit(1);
        }

        ip = (int *)calloc((size_t)(2 * sqrt(np) + 4), sizeof(int));
        if (ip == NULL) {
            printf("We could not obtain enough memory for ip.\n");
            exit(1);
        }

        if ((psvorsh == 3) || (psvorsh == 1)) {
            ip[0] = 0;
            cdft(2 * np, 1, psvt, ip, w);
        }
        if ((psvorsh == 3) || (psvorsh == 2)) {
            ip[0] = 0;
            cdft(2 * np, 1, sht, ip, w);
        }

        free(ip);
        free(w);

        x = (double *)calloc((size_t)np, sizeof(double));
        if (x == NULL) {
            printf("We could not obtain enough memory for x.\n");
            exit(1);
        }

        ypsvt = (double *)calloc((size_t)np, sizeof(double));
        if (ypsvt == NULL) {
            printf("We could not obtain enough memory for ypsvt.\n");
            exit(1);
        }

        ysht = (double *)calloc((size_t)np, sizeof(double));
        if (ysht == NULL) {
            printf("We could not obtain enough memory for ysht.\n");
            exit(1);
        }

        ynett = (double *)calloc((size_t)np, sizeof(double));
        if (ynett == NULL) {
            printf("We could not obtain enough memory for ynett.\n");
            exit(1);
        }

        for (i = 0; i < np; i++) {
            x[i] = (double)(tlen * (double)i / (double)np);
            if (psvorsh != 2)
                ypsvt[i] = (double)(psvt[2 * i] * exp(omegai * (double)x[i]) / tlen * 1000.0);
            if (psvorsh != 1)
                ysht[i] = (double)(sht[2 * i] * exp(omegai * (double)x[i]) / tlen * 1000.0);
            ynett[i] = ypsvt[i] + ysht[i];
        }

        if (psvorsh != 2) free(psvt);
        if (psvorsh != 1) free(sht);

        newsacheader(&hv);
        hv.NPTS = np;
        hv.B = x[0];
        hv.E = x[np - 1];
        hv.DELTA = (double)tlen / np;
        hv.STLA = stla;
        hv.STLO = stlo;
        hv.EVLA = evla;
        hv.EVLO = evlo;
        hv.EVDP = r0;
        strcpy(hv.KCMPNM, "radial");

        if (psvorsh == 3) {
            sprintf(output, "%s%s", spcfile, backr);
            if (asciioutput == 0) {
                writesac(output, &hv, ynett);
            } else if (asciioutput == 1) {
                writeasciisac(output, &hv, ynett);
            }
        }
        if (psvorsh == 1) {
            sprintf(output, "%s%s", spcfile, backr);
            if (asciioutput == 0) {
                writesac(output, &hv, ypsvt);
            } else if (asciioutput == 1) {
                writeasciisac(output, &hv, ypsvt);
            }
        }
        if (psvorsh == 2) {
            sprintf(output, "%s%s", spcfile, backr);
            if (asciioutput == 0) {
                writesac(output, &hv, ysht);
            } else if (asciioutput == 1) {
                writeasciisac(output, &hv, ysht);
            }
        }

        free(x);
        free(ypsvt);
        free(ysht);
        free(ynett);
    } else {
        if (psvorsh != 2) free(psvt);
        if (psvorsh != 1) free(sht);
    }

    if ((component == 7) || (component == 6) || (component == 5) || (component == 4)) {
        w = (double *)calloc((size_t)np * 5 / 2, sizeof(double));
        if (w == NULL) {
            printf("We could not obtain enough memory for w.\n");
            exit(1);
        }

        ip = (int *)calloc((size_t)(2 * sqrt(np) + 4), sizeof(int));
        if (ip == NULL) {
            printf("We could not obtain enough memory for ip.\n");
            exit(1);
        }

        if ((psvorsh == 3) || (psvorsh == 1)) {
            ip[0] = 0;
            cdft(2 * np, 1, psvp, ip, w);
        }
        if ((psvorsh == 3) || (psvorsh == 2)) {
            ip[0] = 0;
            cdft(2 * np, 1, shp, ip, w);
        }

        free(ip);
        free(w);

        x = (double *)calloc((size_t)np, sizeof(double));
        if (x == NULL) {
            printf("We could not obtain enough memory for x.\n");
            exit(1);
        }

        ypsvp = (double *)calloc((size_t)np, sizeof(double));
        if (ypsvp == NULL) {
            printf("We could not obtain enough memory for ypsvp.\n");
            exit(1);
        }

        yshp = (double *)calloc((size_t)np, sizeof(double));
        if (yshp == NULL) {
            printf("We could not obtain enough memory for yshp.\n");
            exit(1);
        }

        ynetp = (double *)calloc((size_t)np, sizeof(double));
        if (ynetp == NULL) {
            printf("We could not obtain enough memory for ynetp.\n");
            exit(1);
        }

        for (i = 0; i < np; i++) {
            x[i] = (double)(tlen * (double)i / (double)np);
            if (psvorsh != 2)
                ypsvp[i] = (double)(psvp[2 * i] * exp(omegai * (double)x[i]) / tlen * 1000.0);
            if (psvorsh != 1)
                yshp[i] = (double)(shp[2 * i] * exp(omegai * (double)x[i]) / tlen * 1000.0);
            ynetp[i] = ypsvp[i] + yshp[i];
        }

        if (psvorsh != 2) free(psvp);
        if (psvorsh != 1) free(shp);

        newsacheader(&hv);
        hv.NPTS = np;
        hv.B = x[0];
        hv.E = x[np - 1];
        hv.DELTA = (double)tlen / np;
        hv.STLA = stla;
        hv.STLO = stlo;
        hv.EVLA = evla;
        hv.EVLO = evlo;
        hv.EVDP = r0;
        strcpy(hv.KCMPNM, "trnsvers");

        if (psvorsh == 3) {
            sprintf(output, "%s%s", spcfile, backt);
            if (asciioutput == 0) {
                writesac(output, &hv, ynetp);
            } else if (asciioutput == 1) {
                writeasciisac(output, &hv, ynetp);
            }
        }
        if (psvorsh == 1) {
            sprintf(output, "%s%s", spcfile, backt);
            if (asciioutput == 0) {
                writesac(output, &hv, ypsvp);
            } else if (asciioutput == 1) {
                writeasciisac(output, &hv, ypsvp);
            }
        }
        if (psvorsh == 2) {
            sprintf(output, "%s%s", spcfile, backt);
            if (asciioutput == 0) {
                writesac(output, &hv, yshp);
            } else if (asciioutput == 1) {
                writeasciisac(output, &hv, yshp);
            }
        }

        free(x);
        free(ypsvp);
        free(yshp);
        free(ynetp);
    } else {
        if (psvorsh != 2) free(psvp);
        if (psvorsh != 1) free(shp);
    }
}

void gcc2ggc(double *theta) {
    const double fl = 0.00335293;
    double ec2;
    ec2 = 2.0 * fl - fl * fl;
    double onemec2 = 1.0 - ec2;
    *theta = 90.0 - atan(tan((90.0 - *theta) / 180.0 * PI) / onemec2) * 180.0 / PI;
}

void helpmessage(void) {
    puts("USAGE:  spcsac [option]");
    puts("-v:            See the version of this software.");
    puts("-a:            ASCII output for non-SAC-users.");
    puts("-h:            See the help message (this text).");
    puts("-e:            If you want to work in a prompt mode rather than");
    puts("               in a command line mode.");
    puts("-l LSMOOTH:    Set the smoothness parameter to the specified value LSMOOTH.");
    puts("               If you want to know the appropriate lsmooth for your");
    puts("               expected sampling frequency,");
    puts("               you can calculate with this spcsac software:");
    puts("               %spcsac -i");
    puts("-i:            In order to obtain appropriate lsmooth, spcsac will");
    puts("               help you to calculate the value.");
    puts("-H SAMPLINGHZ: If you put Sampling Hz value, spcsac will automatically");
    puts("               calculate the appropriate LSMOOTH value.");
    puts("-d DIRECTORY:  All spc files in DIRECTORY will be transformed.");
    puts("               default: (./)");
    puts("-f FILE:       One specified spc FILE only will be transformed.");
    puts("-p/-s STRINGS: Set the suffixes that the input SH and PSV spc files shall have.");
    puts("               (default: .PSV.spc and .SH.spc");
    puts("-c RTZ:        You can choose which component to be transformed.");
    puts("               (R/T/Z)               (default: RTZ)");
    puts("-m MODE:       You can choose PSV/SH channel. (SH/PSV, sh/psv) (default: PSVSH)");
    puts("-r/-t/-z STR:  Set the suffixes of the radial/transverse/vertical");
    puts("               component output SAC files. (default: .Rs/.Ts/.Zs)");
    puts("");
    puts("If you do not put any options, you are to transform all sets of");
    puts(" PSV and SH spc files named \"xxxxxx.PSV.spc\" and \"xxxxxx.SH.spc\"");
    puts(" in the current working directory. LSMOOTH is set as 4.");
    puts(" You will obtain complete seismograms with names like \"xxxxxx.Rs\".");
    puts("");
    puts("EXAMPLES:");
    puts("\tspcsac -l 8 -mPSV -ppsv.spc_tmp -cZT -t .transverse.sac");
    puts("");
    puts("Check the spcsac manual for more information.");
    exit(1);
}

int filexist(const char *filename) {
    FILE *fp;

    if ((fp = fopen(filename, "r")) == NULL) return (0);
    fclose(fp);
    return (1);
}

void lsmoothfinder(void) {
    float tlen;
    int np, np0;
    float freq;
    int lsmooth, i;
    FILE *fp;
    char tmpchar[80];
    puts("This is LSMOOTH finder.\n");
    strcpy(tmpchar, "\0");
    printf("Which scp file?");
    scanf("%s", tmpchar);
    if (!(fp = fopen(tmpchar, "r"))) error01(tmpchar);
    if (!fscanf(fp, "%f\n%d\n", &tlen, &np0)) error01(tmpchar);
    fclose(fp);
    np = 1;
    while (np < np0) np = np * 2;
    printf("How much do you expect as the sampling frequency? (in Hz)\t");
    scanf("%f", &freq);
    lsmooth = (int)(0.5 * tlen * freq / (float)np);
    i = 1;
    while (i < lsmooth) i = i * 2;
    lsmooth = i;
    printf("You should set lsmooth as %d\n", lsmooth);
    printf("You can realize this by putting -l option.\n");
    printf("\te.g. %s spcsac -l %d\n", "%", lsmooth);
    exit(1);
}

int lsmoothfinder_(float tlen, int np0, float freq) {
    int np;
    int lsmooth, i;
    FILE *fp;
    char tmpchar[80];
    np = 1;
    while (np < np0) np = np * 2;
    lsmooth = (int)(0.5 * tlen * freq / (float)np);
    i = 1;
    while (i < lsmooth) i = i * 2;
    lsmooth = i;
    return lsmooth;
}

void nojokenoscience(void) {
    int mode;
    srand((unsigned)time(NULL));
    mode = (rand() % 18 + 1);
    if (mode == 1) {
        puts("NO CAFFEINE, NO SCIENCE!");
    } else if (mode == 2) {
        puts("A girl at a loss asked a policeman in New York,");
        puts("  \"how can I get to Carnegie Hall?\"");
        puts("  \"Practice,\" answered the policeman.");
    } else if (mode == 3) {
        puts("People don't like NIH.");
        puts("......................");
        puts("Not Invented Here.");
    } else if (mode == 4) {
        puts("We have almost completely established earthquake prediction method.");
        puts("  We can predict 100% of earthquakes which took place two seconds ago.");
        puts("  The only thing we have to do is now pure and simple:");
        puts("  We have to advance the timing to predict four more seconds.");
    } else if (mode == 5) {
        puts("Brazil is the country of the future:");
        puts("       and forever it will be.");
    } else if (mode == 6) {
        puts("\"What is your opinion about the shortage of meats?\"\n");
        puts("   \"What do you mean by \'meats\'?,\" replied a North-Korean,");
        puts("   \"I can't figure out what \'your opinion\' means,\" wondered a Chinese,");
        puts("   \"I've never heard of the word \'shortage\' before!\" screamed an American.");
    } else if (mode == 7) {
        puts("From IBM Fortran Compiler Manual(1977)");
        puts("Error Code 103: Correct error and resubmit your problem.");
    } else if (mode == 8) {
        puts("Vice President Richard B. Cheney:");
        puts("He calls the shots.\n");
        puts("\"The message is, if Dick Cheney is willing to shoot an innocent American");
        puts("   citizen at point-blank range, imagine what he'll do to you,\"");
        puts("Mr. Bush said.");
    } else if (mode == 9) {
        puts("Trust me but verify!");
    } else if (mode == 10) {
        puts("Poacher poached.");
    } else if (mode == 11) {
        puts("There's something we hadn't better know.");
    } else if (mode == 12) {
        puts("A man felt around a lamp stand for his car key one night.");
        puts("His friend came up and asked:");
        puts("   \"What are you doing? \"");
        puts("   \"I've lost my car key and I'm looking for it,\" answered the man.");
        puts("   \"You've lost it around here?");
        puts("   \"No. Up there.\"");
        puts("   \"Then how come you so eager about searching here? \"");
        puts("   \"Look. Here it's brighter and easier to search.\"");
    } else if (mode == 13) {
        puts("And so..... what's the xxxxing result?");
    } else if (mode == 14) {
        puts("Surely You're Joking, Mr. Geller!");
    } else if (mode == 15) {
        puts("There are three kinds of people in the world:");
        puts("   those who can count and those who cannot.");
    } else if (mode == 16) {
        puts("No doubt Paul Dirac was a great physicist.");
        puts("But of course he could make a mistake.");
        puts("In a class of physics his student asked Dr. Dirac,");
        puts("    \"Sir, I think that you had a mistake in writing,\"");
        puts("    \"That \'c\' squared must be \'c\' itself, I think.\"");
        puts("Prof. Dirac glanced at the blackboad and replied,");
        puts("    \"Well, what is your definition of \'c\'?");
        puts("    \"Velocity of light, sir,\" replied the student immediately.");
        puts("    \"Oh, I'm sorry but my definition of the velocity of light is.....\"");
        puts("    \"\'c\' squared.\"");
        puts("He continued the class for two more hours.");
        puts("");
        puts(" with his \'original\' definition of the velocity of light.");
    } else if (mode == 17) {
        puts("We will rebuilt the detector, there's no question.");
        puts("                                   --- Yoji Totsuka");
    } else {
        puts("You have finally found me! ");
        puts("Thank you very much for working with DSM softwares.");
        puts("DSM was developed by the Global Seismology Group at the Univ. of Tokyo.");
        puts("You can see our site here:");
        puts("https://utglobalseismology.org");
    }
    exit(1);
}

void error01(char *spcfile) {
    printf("spcfile %s is not appropriate.\nPlease check it and if you want to know more,\n", spcfile);
    puts("you can read the help message by typing");
    puts("> spcsac -h");
    exit(1);
}

void modeprompt(int *filordir, int *psvorsh, int *component, char *backr,
                char *backt, char *backz, char *backp, char *backs,
                char *spcfile) {
    char tmpchar[80];
    strcpy(tmpchar, "");
    puts("-----------------------------------------------------------------------");
    puts("This is spcsac command prompt mode. You are to be asked some questions:");
    puts("When encountering the question which you do not care, which you do not");
    puts("even have to think of or whose default answer you agree,");
    puts("please press \"\\\" key (backslash) in order to skip the question.");
    puts("-----------------------------------------------------------------------");
    puts("Do you want to transform an spc file");
    puts("\trather than all spc files in a certain directory?");
    printf("(y/n/q) (default n)");
    scanf("%s", tmpchar);
    if (!strcmp(tmpchar, "\\")) {
        strcpy(tmpchar, "n");
    }
    if (!strcmp(tmpchar, "y")) {
        *filordir = 1;
        printf("Which file do you want to transform?\t");
        scanf("%s", tmpchar);
        strcpy(spcfile, tmpchar);
    } else if (!strcmp(tmpchar, "n")) {
        printf("What directory is your target device? (default: ./)");
        scanf("%s", tmpchar);
        if (strcmp(tmpchar, "\\")) strcpy(spcfile, tmpchar);
        puts("Do you want to obtain complete (PSV+SH) seismograms or not?");
        printf("(y/n) (default y)");
        strcpy(tmpchar, "y");
        scanf("%s", tmpchar);
        if (!strcmp(tmpchar, "n")) {
            printf("Which seismogram? (PSV/SH)");
            scanf("%s", tmpchar);
            if (!strcmp(tmpchar, "PSV")) {
                *psvorsh = 1;
            } else if (!strcmp(tmpchar, "SH")) {
                *psvorsh = 2;
            } else {
                puts("Your answer was not appropriate, but we are going to make complete ones.");
            }
        } else if (!strcmp(tmpchar, "y")) {
        } else {
            puts("Your answer was not appropriate, but we are going to make complete ones.");
        }

        if ((*psvorsh == 3) || (*psvorsh == 1)) {
            puts("How do you distinguish PSV spc file? (default: .PSV.spc)");
            printf("You are to required to put common strings at the back of file names.");
            scanf("%s", tmpchar);
            if (strcmp(tmpchar, "\\")) strcpy(backp, tmpchar);
        }
        if ((*psvorsh == 3) || (*psvorsh == 2)) {
            puts("How do you distinguish SH spc file? (default: .SH.spc)");
            printf("You are to required to put common strings at the back of file names.");
            scanf("%s", tmpchar);
            if (strcmp(tmpchar, "\\")) strcpy(backs, tmpchar);
        }
    } else if (!strcmp(tmpchar, "q")) {
        exit(1);
    } else {
        puts("spcsac command prompt mode error occurred.");
        helpmessage();
    }
    puts("Which component of seismograms do you want to obtain? Select from R/T/Z.");
    printf("\t(default RTZ)");
    scanf("%s", tmpchar);
    if (!strcmp(tmpchar, "T")) *component = 4;
    if (!strcmp(tmpchar, "R")) *component = 2;
    if (!strcmp(tmpchar, "Z")) *component = 1;
    if (!strcmp(tmpchar, "TR")) *component = 6;
    if (!strcmp(tmpchar, "RT")) *component = 6;
    if (!strcmp(tmpchar, "TZ")) *component = 5;
    if (!strcmp(tmpchar, "ZT")) *component = 5;
    if (!strcmp(tmpchar, "RZ")) *component = 3;
    if (!strcmp(tmpchar, "ZR")) *component = 3;
    if (*filordir == 0) {
        if ((*component == 7) || (*component == 6) || (*component == 3) ||
            (*component == 2)) {
            puts("How do you distinguish R sac files? (default: .Rs)");
            printf("You are to required to put common strings at the back of file names.");
            scanf("%s", tmpchar);
            if (strcmp(tmpchar, "\\")) strcpy(backr, tmpchar);
        }
        if ((*component == 7) || (*component == 6) || (*component == 5) ||
            (*component == 4)) {
            puts("How do you distinguish T sac files? (default: .Ts)");
            printf("You are to required to put common strings at the back of file names.");
            scanf("%s", tmpchar);
            if (strcmp(tmpchar, "\\")) strcpy(backt, tmpchar);
        }
        if ((*component == 7) || (*component == 5) || (*component == 3) ||
            (*component == 1)) {
            puts("How do you distinguish Z sac files? (default: .Zs)");
            printf("You are to required to put commmon strings at the back of file names.");
            scanf("%s", tmpchar);
            if (strcmp(tmpchar, "\\")) strcpy(backz, tmpchar);
        }
    }

    if (*filordir == 1) {
        if ((*component == 7) || (*component == 6) || (*component == 3) ||
            (*component == 2)) {
            printf("What is the name of radial component sac file? (default: tmpsac.Rs)");
            scanf("%s", tmpchar);
            if (strcmp(tmpchar, "\\")) strcpy(backr, tmpchar);
            if (!strcmp(tmpchar, "\\")) strcpy(backr, "tmpsac.Rs");
        }
        if ((*component == 7) || (*component == 6) || (*component == 5) ||
            (*component == 4)) {
            printf("What is the name of transverse component sac file? (default: tmpsac.Ts)");
            scanf("%s", tmpchar);
            if (strcmp(tmpchar, "\\")) strcpy(backt, tmpchar);
            if (!strcmp(tmpchar, "\\")) strcpy(backt, "tmpsac.Ts");
        }
        if ((*component == 7) || (*component == 5) || (*component == 3) ||
            (*component == 1)) {
            printf("What is the name of vertical component sac file? (default: tmpsac.Zs)");
            scanf("%s", tmpchar);
            if (strcmp(tmpchar, "\\")) strcpy(backz, tmpchar);
            if (!strcmp(tmpchar, "\\")) strcpy(backz, "tmpsac.Zs");
        }
    }
}

void cdft(long n, int isgn, double *a, int *ip, double *w) {
    void makewt(long nw, int *ip, double *w);
    void bitrv2(long n, int *ip, double *a);
    void bitrv2conj(long n, int *ip, double *a);
    void cftfsub(long n, double *a, double *w);
    void cftbsub(long n, double *a, double *w);

    if (n > (ip[0] << 2)) {
        makewt(n >> 2, ip, w);
    }
    if (n > 4) {
        if (isgn >= 0) {
            bitrv2(n, ip + 2, a);
            cftfsub(n, a, w);
        } else {
            bitrv2conj(n, ip + 2, a);
            cftbsub(n, a, w);
        }
    } else if (n == 4) {
        cftfsub(n, a, w);
    }
}

void rdft(long n, int isgn, double *a, int *ip, double *w) {
    void makewt(long nw, int *ip, double *w);
    void makect(long nc, int *ip, double *c);
    void bitrv2(long n, int *ip, double *a);
    void cftfsub(long n, double *a, double *w);
    void cftbsub(long n, double *a, double *w);
    void rftfsub(long n, double *a, int nc, double *c);
    void rftbsub(long n, double *a, int nc, double *c);
    long nw, nc;
    double xi;

    nw = ip[0];
    if (n > (nw << 2)) {
        nw = n >> 2;
        makewt(nw, ip, w);
    }
    nc = ip[1];
    if (n > (nc << 2)) {
        nc = n >> 2;
        makect(nc, ip, w + nw);
    }
    if (isgn >= 0) {
        if (n > 4) {
            bitrv2(n, ip + 2, a);
            cftfsub(n, a, w);
            rftfsub(n, a, nc, w + nw);
        } else if (n == 4) {
            cftfsub(n, a, w);
        }
        xi = a[0] - a[1];
        a[0] += a[1];
        a[1] = xi;
    } else {
        a[1] = 0.5 * (a[0] - a[1]);
        a[0] -= a[1];
        if (n > 4) {
            rftbsub(n, a, nc, w + nw);
            bitrv2(n, ip + 2, a);
            cftbsub(n, a, w);
        } else if (n == 4) {
            cftfsub(n, a, w);
        }
    }
}

void ddct(long n, int isgn, double *a, int *ip, double *w) {
    void makewt(long nw, int *ip, double *w);
    void makect(long nc, int *ip, double *c);
    void bitrv2(long n, int *ip, double *a);
    void cftfsub(long n, double *a, double *w);
    void cftbsub(long n, double *a, double *w);
    void rftfsub(long n, double *a, int nc, double *c);
    void rftbsub(long n, double *a, int nc, double *c);
    void dctsub(long n, double *a, int nc, double *c);
    long j, nw, nc;
    double xr;

    nw = ip[0];
    if (n > (nw << 2)) {
        nw = n >> 2;
        makewt(nw, ip, w);
    }
    nc = ip[1];
    if (n > nc) {
        nc = n;
        makect(nc, ip, w + nw);
    }
    if (isgn < 0) {
        xr = a[n - 1];
        for (j = n - 2; j >= 2; j -= 2) {
            a[j + 1] = a[j] - a[j - 1];
            a[j] += a[j - 1];
        }
        a[1] = a[0] - xr;
        a[0] += xr;
        if (n > 4) {
            rftbsub(n, a, nc, w + nw);
            bitrv2(n, ip + 2, a);
            cftbsub(n, a, w);
        } else if (n == 4) {
            cftfsub(n, a, w);
        }
    }
    dctsub(n, a, nc, w + nw);
    if (isgn >= 0) {
        if (n > 4) {
            bitrv2(n, ip + 2, a);
            cftfsub(n, a, w);
            rftfsub(n, a, nc, w + nw);
        } else if (n == 4) {
            cftfsub(n, a, w);
        }
        xr = a[0] - a[1];
        a[0] += a[1];
        for (j = 2; j < n; j += 2) {
            a[j - 1] = a[j] - a[j + 1];
            a[j] += a[j + 1];
        }
        a[n - 1] = xr;
    }
}

void ddst(long n, int isgn, double *a, int *ip, double *w) {
    void makewt(long nw, int *ip, double *w);
    void makect(long nc, int *ip, double *c);
    void bitrv2(long n, int *ip, double *a);
    void cftfsub(long n, double *a, double *w);
    void cftbsub(long n, double *a, double *w);
    void rftfsub(long n, double *a, int nc, double *c);
    void rftbsub(long n, double *a, int nc, double *c);
    void dstsub(long n, double *a, int nc, double *c);
    long j, nw, nc;
    double xr;

    nw = ip[0];
    if (n > (nw << 2)) {
        nw = n >> 2;
        makewt(nw, ip, w);
    }
    nc = ip[1];
    if (n > nc) {
        nc = n;
        makect(nc, ip, w + nw);
    }
    if (isgn < 0) {
        xr = a[n - 1];
        for (j = n - 2; j >= 2; j -= 2) {
            a[j + 1] = -a[j] - a[j - 1];
            a[j] -= a[j - 1];
        }
        a[1] = a[0] + xr;
        a[0] -= xr;
        if (n > 4) {
            rftbsub(n, a, nc, w + nw);
            bitrv2(n, ip + 2, a);
            cftbsub(n, a, w);
        } else if (n == 4) {
            cftfsub(n, a, w);
        }
    }
    dstsub(n, a, nc, w + nw);
    if (isgn >= 0) {
        if (n > 4) {
            bitrv2(n, ip + 2, a);
            cftfsub(n, a, w);
            rftfsub(n, a, nc, w + nw);
        } else if (n == 4) {
            cftfsub(n, a, w);
        }
        xr = a[0] - a[1];
        a[0] += a[1];
        for (j = 2; j < n; j += 2) {
            a[j - 1] = -a[j] - a[j + 1];
            a[j] -= a[j + 1];
        }
        a[n - 1] = -xr;
    }
}

void dfct(long n, double *a, double *t, int *ip, double *w) {
    void makewt(long nw, int *ip, double *w);
    void makect(long nc, int *ip, double *c);
    void bitrv2(long n, int *ip, double *a);
    void cftfsub(long n, double *a, double *w);
    void rftfsub(long n, double *a, int nc, double *c);
    void dctsub(long n, double *a, int nc, double *c);
    long j, k, l, m, mh, nw, nc;
    double xr, xi, yr, yi;

    nw = ip[0];
    if (n > (nw << 3)) {
        nw = n >> 3;
        makewt(nw, ip, w);
    }
    nc = ip[1];
    if (n > (nc << 1)) {
        nc = n >> 1;
        makect(nc, ip, w + nw);
    }
    m = n >> 1;
    yi = a[m];
    xi = a[0] + a[n];
    a[0] -= a[n];
    t[0] = xi - yi;
    t[m] = xi + yi;
    if (n > 2) {
        mh = m >> 1;
        for (j = 1; j < mh; j++) {
            k = m - j;
            xr = a[j] - a[n - j];
            xi = a[j] + a[n - j];
            yr = a[k] - a[n - k];
            yi = a[k] + a[n - k];
            a[j] = xr;
            a[k] = yr;
            t[j] = xi - yi;
            t[k] = xi + yi;
        }
        t[mh] = a[mh] + a[n - mh];
        a[mh] -= a[n - mh];
        dctsub(m, a, nc, w + nw);
        if (m > 4) {
            bitrv2(m, ip + 2, a);
            cftfsub(m, a, w);
            rftfsub(m, a, nc, w + nw);
        } else if (m == 4) {
            cftfsub(m, a, w);
        }
        a[n - 1] = a[0] - a[1];
        a[1] = a[0] + a[1];
        for (j = m - 2; j >= 2; j -= 2) {
            a[2 * j + 1] = a[j] + a[j + 1];
            a[2 * j - 1] = a[j] - a[j + 1];
        }
        l = 2;
        m = mh;
        while (m >= 2) {
            dctsub(m, t, nc, w + nw);
            if (m > 4) {
                bitrv2(m, ip + 2, t);
                cftfsub(m, t, w);
                rftfsub(m, t, nc, w + nw);
            } else if (m == 4) {
                cftfsub(m, t, w);
            }
            a[n - l] = t[0] - t[1];
            a[l] = t[0] + t[1];
            k = 0;
            for (j = 2; j < m; j += 2) {
                k += l << 2;
                a[k - l] = t[j] - t[j + 1];
                a[k + l] = t[j] + t[j + 1];
            }
            l <<= 1;
            mh = m >> 1;
            for (j = 0; j < mh; j++) {
                k = m - j;
                t[j] = t[m + k] - t[m + j];
                t[k] = t[m + k] + t[m + j];
            }
            t[mh] = t[m + mh];
            m = mh;
        }
        a[l] = t[0];
        a[n] = t[2] - t[1];
        a[0] = t[2] + t[1];
    } else {
        a[1] = a[0];
        a[2] = t[0];
        a[0] = t[1];
    }
}

void dfst(long n, double *a, double *t, int *ip, double *w) {
    void makewt(long nw, int *ip, double *w);
    void makect(long nc, int *ip, double *c);
    void bitrv2(long n, int *ip, double *a);
    void cftfsub(long n, double *a, double *w);
    void rftfsub(long n, double *a, int nc, double *c);
    void dstsub(long n, double *a, int nc, double *c);
    long j, k, l, m, mh, nw, nc;
    double xr, xi, yr, yi;

    nw = ip[0];
    if (n > (nw << 3)) {
        nw = n >> 3;
        makewt(nw, ip, w);
    }
    nc = ip[1];
    if (n > (nc << 1)) {
        nc = n >> 1;
        makect(nc, ip, w + nw);
    }
    if (n > 2) {
        m = n >> 1;
        mh = m >> 1;
        for (j = 1; j < mh; j++) {
            k = m - j;
            xr = a[j] + a[n - j];
            xi = a[j] - a[n - j];
            yr = a[k] + a[n - k];
            yi = a[k] - a[n - k];
            a[j] = xr;
            a[k] = yr;
            t[j] = xi + yi;
            t[k] = xi - yi;
        }
        t[0] = a[mh] - a[n - mh];
        a[mh] += a[n - mh];
        a[0] = a[m];
        dstsub(m, a, nc, w + nw);
        if (m > 4) {
            bitrv2(m, ip + 2, a);
            cftfsub(m, a, w);
            rftfsub(m, a, nc, w + nw);
        } else if (m == 4) {
            cftfsub(m, a, w);
        }
        a[n - 1] = a[1] - a[0];
        a[1] = a[0] + a[1];
        for (j = m - 2; j >= 2; j -= 2) {
            a[2 * j + 1] = a[j] - a[j + 1];
            a[2 * j - 1] = -a[j] - a[j + 1];
        }
        l = 2;
        m = mh;
        while (m >= 2) {
            dstsub(m, t, nc, w + nw);
            if (m > 4) {
                bitrv2(m, ip + 2, t);
                cftfsub(m, t, w);
                rftfsub(m, t, nc, w + nw);
            } else if (m == 4) {
                cftfsub(m, t, w);
            }
            a[n - l] = t[1] - t[0];
            a[l] = t[0] + t[1];
            k = 0;
            for (j = 2; j < m; j += 2) {
                k += l << 2;
                a[k - l] = -t[j] - t[j + 1];
                a[k + l] = t[j] - t[j + 1];
            }
            l <<= 1;
            mh = m >> 1;
            for (j = 1; j < mh; j++) {
                k = m - j;
                t[j] = t[m + k] + t[m + j];
                t[k] = t[m + k] - t[m + j];
            }
            t[0] = t[m + mh];
            m = mh;
        }
        a[l] = t[0];
    }
    a[0] = 0;
}

/* -------- initializing routines -------- */

#include <math.h>

void makewt(long nw, int *ip, double *w) {
    void bitrv2(long n, int *ip, double *a);
    long j, nwh;
    double delta, x, y;

    ip[0] = nw;
    ip[1] = 1;
    if (nw > 2) {
        nwh = nw >> 1;
        delta = atan(1.0) / nwh;
        w[0] = 1;
        w[1] = 0;
        w[nwh] = cos(delta * nwh);
        w[nwh + 1] = w[nwh];
        if (nwh > 2) {
            for (j = 2; j < nwh; j += 2) {
                x = cos(delta * j);
                y = sin(delta * j);
                w[j] = x;
                w[j + 1] = y;
                w[nw - j] = y;
                w[nw - j + 1] = x;
            }
            bitrv2(nw, ip + 2, w);
        }
    }
}

void makect(long nc, int *ip, double *c) {
    long j, nch;
    double delta;

    ip[1] = nc;
    if (nc > 1) {
        nch = nc >> 1;
        delta = atan(1.0) / nch;
        c[0] = cos(delta * nch);
        c[nch] = 0.5 * c[0];
        for (j = 1; j < nch; j++) {
            c[j] = 0.5 * cos(delta * j);
            c[nc - j] = 0.5 * sin(delta * j);
        }
    }
}

/* -------- child routines -------- */

void bitrv2(long n, int *ip, double *a) {
    long j, j1, k, k1, l, m, m2;
    double xr, xi, yr, yi;

    ip[0] = 0;
    l = n;
    m = 1;
    while ((m << 3) < l) {
        l >>= 1;
        for (j = 0; j < m; j++) {
            ip[m + j] = ip[j] + l;
        }
        m <<= 1;
    }
    m2 = 2 * m;
    if ((m << 3) == l) {
        for (k = 0; k < m; k++) {
            for (j = 0; j < k; j++) {
                j1 = 2 * j + ip[k];
                k1 = 2 * k + ip[j];
                xr = a[j1];
                xi = a[j1 + 1];
                yr = a[k1];
                yi = a[k1 + 1];
                a[j1] = yr;
                a[j1 + 1] = yi;
                a[k1] = xr;
                a[k1 + 1] = xi;
                j1 += m2;
                k1 += 2 * m2;
                xr = a[j1];
                xi = a[j1 + 1];
                yr = a[k1];
                yi = a[k1 + 1];
                a[j1] = yr;
                a[j1 + 1] = yi;
                a[k1] = xr;
                a[k1 + 1] = xi;
                j1 += m2;
                k1 -= m2;
                xr = a[j1];
                xi = a[j1 + 1];
                yr = a[k1];
                yi = a[k1 + 1];
                a[j1] = yr;
                a[j1 + 1] = yi;
                a[k1] = xr;
                a[k1 + 1] = xi;
                j1 += m2;
                k1 += 2 * m2;
                xr = a[j1];
                xi = a[j1 + 1];
                yr = a[k1];
                yi = a[k1 + 1];
                a[j1] = yr;
                a[j1 + 1] = yi;
                a[k1] = xr;
                a[k1 + 1] = xi;
            }
            j1 = 2 * k + m2 + ip[k];
            k1 = j1 + m2;
            xr = a[j1];
            xi = a[j1 + 1];
            yr = a[k1];
            yi = a[k1 + 1];
            a[j1] = yr;
            a[j1 + 1] = yi;
            a[k1] = xr;
            a[k1 + 1] = xi;
        }
    } else {
        for (k = 1; k < m; k++) {
            for (j = 0; j < k; j++) {
                j1 = 2 * j + ip[k];
                k1 = 2 * k + ip[j];
                xr = a[j1];
                xi = a[j1 + 1];
                yr = a[k1];
                yi = a[k1 + 1];
                a[j1] = yr;
                a[j1 + 1] = yi;
                a[k1] = xr;
                a[k1 + 1] = xi;
                j1 += m2;
                k1 += m2;
                xr = a[j1];
                xi = a[j1 + 1];
                yr = a[k1];
                yi = a[k1 + 1];
                a[j1] = yr;
                a[j1 + 1] = yi;
                a[k1] = xr;
                a[k1 + 1] = xi;
            }
        }
    }
}

void bitrv2conj(long n, int *ip, double *a) {
    long j, j1, k, k1, l, m, m2;
    double xr, xi, yr, yi;

    ip[0] = 0;
    l = n;
    m = 1;
    while ((m << 3) < l) {
        l >>= 1;
        for (j = 0; j < m; j++) {
            ip[m + j] = ip[j] + l;
        }
        m <<= 1;
    }
    m2 = 2 * m;
    if ((m << 3) == l) {
        for (k = 0; k < m; k++) {
            for (j = 0; j < k; j++) {
                j1 = 2 * j + ip[k];
                k1 = 2 * k + ip[j];
                xr = a[j1];
                xi = -a[j1 + 1];
                yr = a[k1];
                yi = -a[k1 + 1];
                a[j1] = yr;
                a[j1 + 1] = yi;
                a[k1] = xr;
                a[k1 + 1] = xi;
                j1 += m2;
                k1 += 2 * m2;
                xr = a[j1];
                xi = -a[j1 + 1];
                yr = a[k1];
                yi = -a[k1 + 1];
                a[j1] = yr;
                a[j1 + 1] = yi;
                a[k1] = xr;
                a[k1 + 1] = xi;
                j1 += m2;
                k1 -= m2;
                xr = a[j1];
                xi = -a[j1 + 1];
                yr = a[k1];
                yi = -a[k1 + 1];
                a[j1] = yr;
                a[j1 + 1] = yi;
                a[k1] = xr;
                a[k1 + 1] = xi;
                j1 += m2;
                k1 += 2 * m2;
                xr = a[j1];
                xi = -a[j1 + 1];
                yr = a[k1];
                yi = -a[k1 + 1];
                a[j1] = yr;
                a[j1 + 1] = yi;
                a[k1] = xr;
                a[k1 + 1] = xi;
            }
            k1 = 2 * k + ip[k];
            a[k1 + 1] = -a[k1 + 1];
            j1 = k1 + m2;
            k1 = j1 + m2;
            xr = a[j1];
            xi = -a[j1 + 1];
            yr = a[k1];
            yi = -a[k1 + 1];
            a[j1] = yr;
            a[j1 + 1] = yi;
            a[k1] = xr;
            a[k1 + 1] = xi;
            k1 += m2;
            a[k1 + 1] = -a[k1 + 1];
        }
    } else {
        a[1] = -a[1];
        a[m2 + 1] = -a[m2 + 1];
        for (k = 1; k < m; k++) {
            for (j = 0; j < k; j++) {
                j1 = 2 * j + ip[k];
                k1 = 2 * k + ip[j];
                xr = a[j1];
                xi = -a[j1 + 1];
                yr = a[k1];
                yi = -a[k1 + 1];
                a[j1] = yr;
                a[j1 + 1] = yi;
                a[k1] = xr;
                a[k1 + 1] = xi;
                j1 += m2;
                k1 += m2;
                xr = a[j1];
                xi = -a[j1 + 1];
                yr = a[k1];
                yi = -a[k1 + 1];
                a[j1] = yr;
                a[j1 + 1] = yi;
                a[k1] = xr;
                a[k1 + 1] = xi;
            }
            k1 = 2 * k + ip[k];
            a[k1 + 1] = -a[k1 + 1];
            a[k1 + m2 + 1] = -a[k1 + m2 + 1];
        }
    }
}

void cftfsub(long n, double *a, double *w) {
    void cft1st(long n, double *a, double *w);
    void cftmdl(long n, int l, double *a, double *w);
    long j, j1, j2, j3, l;
    double x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i;

    l = 2;
    if (n > 8) {
        cft1st(n, a, w);
        l = 8;
        while ((l << 2) < n) {
            cftmdl(n, l, a, w);
            l <<= 2;
        }
    }
    if ((l << 2) == n) {
        for (j = 0; j < l; j += 2) {
            j1 = j + l;
            j2 = j1 + l;
            j3 = j2 + l;
            x0r = a[j] + a[j1];
            x0i = a[j + 1] + a[j1 + 1];
            x1r = a[j] - a[j1];
            x1i = a[j + 1] - a[j1 + 1];
            x2r = a[j2] + a[j3];
            x2i = a[j2 + 1] + a[j3 + 1];
            x3r = a[j2] - a[j3];
            x3i = a[j2 + 1] - a[j3 + 1];
            a[j] = x0r + x2r;
            a[j + 1] = x0i + x2i;
            a[j2] = x0r - x2r;
            a[j2 + 1] = x0i - x2i;
            a[j1] = x1r - x3i;
            a[j1 + 1] = x1i + x3r;
            a[j3] = x1r + x3i;
            a[j3 + 1] = x1i - x3r;
        }
    } else {
        for (j = 0; j < l; j += 2) {
            j1 = j + l;
            x0r = a[j] - a[j1];
            x0i = a[j + 1] - a[j1 + 1];
            a[j] += a[j1];
            a[j + 1] += a[j1 + 1];
            a[j1] = x0r;
            a[j1 + 1] = x0i;
        }
    }
}

void cftbsub(long n, double *a, double *w) {
    void cft1st(long n, double *a, double *w);
    void cftmdl(long n, int l, double *a, double *w);
    long j, j1, j2, j3, l;
    double x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i;

    l = 2;
    if (n > 8) {
        cft1st(n, a, w);
        l = 8;
        while ((l << 2) < n) {
            cftmdl(n, l, a, w);
            l <<= 2;
        }
    }
    if ((l << 2) == n) {
        for (j = 0; j < l; j += 2) {
            j1 = j + l;
            j2 = j1 + l;
            j3 = j2 + l;
            x0r = a[j] + a[j1];
            x0i = -a[j + 1] - a[j1 + 1];
            x1r = a[j] - a[j1];
            x1i = -a[j + 1] + a[j1 + 1];
            x2r = a[j2] + a[j3];
            x2i = a[j2 + 1] + a[j3 + 1];
            x3r = a[j2] - a[j3];
            x3i = a[j2 + 1] - a[j3 + 1];
            a[j] = x0r + x2r;
            a[j + 1] = x0i - x2i;
            a[j2] = x0r - x2r;
            a[j2 + 1] = x0i + x2i;
            a[j1] = x1r - x3i;
            a[j1 + 1] = x1i - x3r;
            a[j3] = x1r + x3i;
            a[j3 + 1] = x1i + x3r;
        }
    } else {
        for (j = 0; j < l; j += 2) {
            j1 = j + l;
            x0r = a[j] - a[j1];
            x0i = -a[j + 1] + a[j1 + 1];
            a[j] += a[j1];
            a[j + 1] = -a[j + 1] - a[j1 + 1];
            a[j1] = x0r;
            a[j1 + 1] = x0i;
        }
    }
}

void cft1st(long n, double *a, double *w) {
    long j, k1, k2;
    double wk1r, wk1i, wk2r, wk2i, wk3r, wk3i;
    double x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i;

    x0r = a[0] + a[2];
    x0i = a[1] + a[3];
    x1r = a[0] - a[2];
    x1i = a[1] - a[3];
    x2r = a[4] + a[6];
    x2i = a[5] + a[7];
    x3r = a[4] - a[6];
    x3i = a[5] - a[7];
    a[0] = x0r + x2r;
    a[1] = x0i + x2i;
    a[4] = x0r - x2r;
    a[5] = x0i - x2i;
    a[2] = x1r - x3i;
    a[3] = x1i + x3r;
    a[6] = x1r + x3i;
    a[7] = x1i - x3r;
    wk1r = w[2];
    x0r = a[8] + a[10];
    x0i = a[9] + a[11];
    x1r = a[8] - a[10];
    x1i = a[9] - a[11];
    x2r = a[12] + a[14];
    x2i = a[13] + a[15];
    x3r = a[12] - a[14];
    x3i = a[13] - a[15];
    a[8] = x0r + x2r;
    a[9] = x0i + x2i;
    a[12] = x2i - x0i;
    a[13] = x0r - x2r;
    x0r = x1r - x3i;
    x0i = x1i + x3r;
    a[10] = wk1r * (x0r - x0i);
    a[11] = wk1r * (x0r + x0i);
    x0r = x3i + x1r;
    x0i = x3r - x1i;
    a[14] = wk1r * (x0i - x0r);
    a[15] = wk1r * (x0i + x0r);
    k1 = 0;
    for (j = 16; j < n; j += 16) {
        k1 += 2;
        k2 = 2 * k1;
        wk2r = w[k1];
        wk2i = w[k1 + 1];
        wk1r = w[k2];
        wk1i = w[k2 + 1];
        wk3r = wk1r - 2 * wk2i * wk1i;
        wk3i = 2 * wk2i * wk1r - wk1i;
        x0r = a[j] + a[j + 2];
        x0i = a[j + 1] + a[j + 3];
        x1r = a[j] - a[j + 2];
        x1i = a[j + 1] - a[j + 3];
        x2r = a[j + 4] + a[j + 6];
        x2i = a[j + 5] + a[j + 7];
        x3r = a[j + 4] - a[j + 6];
        x3i = a[j + 5] - a[j + 7];
        a[j] = x0r + x2r;
        a[j + 1] = x0i + x2i;
        x0r -= x2r;
        x0i -= x2i;
        a[j + 4] = wk2r * x0r - wk2i * x0i;
        a[j + 5] = wk2r * x0i + wk2i * x0r;
        x0r = x1r - x3i;
        x0i = x1i + x3r;
        a[j + 2] = wk1r * x0r - wk1i * x0i;
        a[j + 3] = wk1r * x0i + wk1i * x0r;
        x0r = x1r + x3i;
        x0i = x1i - x3r;
        a[j + 6] = wk3r * x0r - wk3i * x0i;
        a[j + 7] = wk3r * x0i + wk3i * x0r;
        wk1r = w[k2 + 2];
        wk1i = w[k2 + 3];
        wk3r = wk1r - 2 * wk2r * wk1i;
        wk3i = 2 * wk2r * wk1r - wk1i;
        x0r = a[j + 8] + a[j + 10];
        x0i = a[j + 9] + a[j + 11];
        x1r = a[j + 8] - a[j + 10];
        x1i = a[j + 9] - a[j + 11];
        x2r = a[j + 12] + a[j + 14];
        x2i = a[j + 13] + a[j + 15];
        x3r = a[j + 12] - a[j + 14];
        x3i = a[j + 13] - a[j + 15];
        a[j + 8] = x0r + x2r;
        a[j + 9] = x0i + x2i;
        x0r -= x2r;
        x0i -= x2i;
        a[j + 12] = -wk2i * x0r - wk2r * x0i;
        a[j + 13] = -wk2i * x0i + wk2r * x0r;
        x0r = x1r - x3i;
        x0i = x1i + x3r;
        a[j + 10] = wk1r * x0r - wk1i * x0i;
        a[j + 11] = wk1r * x0i + wk1i * x0r;
        x0r = x1r + x3i;
        x0i = x1i - x3r;
        a[j + 14] = wk3r * x0r - wk3i * x0i;
        a[j + 15] = wk3r * x0i + wk3i * x0r;
    }
}

void cftmdl(long n, int l, double *a, double *w) {
    long j, j1, j2, j3, k, k1, k2, m, m2;
    double wk1r, wk1i, wk2r, wk2i, wk3r, wk3i;
    double x0r, x0i, x1r, x1i, x2r, x2i, x3r, x3i;

    m = l << 2;
    for (j = 0; j < l; j += 2) {
        j1 = j + l;
        j2 = j1 + l;
        j3 = j2 + l;
        x0r = a[j] + a[j1];
        x0i = a[j + 1] + a[j1 + 1];
        x1r = a[j] - a[j1];
        x1i = a[j + 1] - a[j1 + 1];
        x2r = a[j2] + a[j3];
        x2i = a[j2 + 1] + a[j3 + 1];
        x3r = a[j2] - a[j3];
        x3i = a[j2 + 1] - a[j3 + 1];
        a[j] = x0r + x2r;
        a[j + 1] = x0i + x2i;
        a[j2] = x0r - x2r;
        a[j2 + 1] = x0i - x2i;
        a[j1] = x1r - x3i;
        a[j1 + 1] = x1i + x3r;
        a[j3] = x1r + x3i;
        a[j3 + 1] = x1i - x3r;
    }
    wk1r = w[2];
    for (j = m; j < l + m; j += 2) {
        j1 = j + l;
        j2 = j1 + l;
        j3 = j2 + l;
        x0r = a[j] + a[j1];
        x0i = a[j + 1] + a[j1 + 1];
        x1r = a[j] - a[j1];
        x1i = a[j + 1] - a[j1 + 1];
        x2r = a[j2] + a[j3];
        x2i = a[j2 + 1] + a[j3 + 1];
        x3r = a[j2] - a[j3];
        x3i = a[j2 + 1] - a[j3 + 1];
        a[j] = x0r + x2r;
        a[j + 1] = x0i + x2i;
        a[j2] = x2i - x0i;
        a[j2 + 1] = x0r - x2r;
        x0r = x1r - x3i;
        x0i = x1i + x3r;
        a[j1] = wk1r * (x0r - x0i);
        a[j1 + 1] = wk1r * (x0r + x0i);
        x0r = x3i + x1r;
        x0i = x3r - x1i;
        a[j3] = wk1r * (x0i - x0r);
        a[j3 + 1] = wk1r * (x0i + x0r);
    }
    k1 = 0;
    m2 = 2 * m;
    for (k = m2; k < n; k += m2) {
        k1 += 2;
        k2 = 2 * k1;
        wk2r = w[k1];
        wk2i = w[k1 + 1];
        wk1r = w[k2];
        wk1i = w[k2 + 1];
        wk3r = wk1r - 2 * wk2i * wk1i;
        wk3i = 2 * wk2i * wk1r - wk1i;
        for (j = k; j < l + k; j += 2) {
            j1 = j + l;
            j2 = j1 + l;
            j3 = j2 + l;
            x0r = a[j] + a[j1];
            x0i = a[j + 1] + a[j1 + 1];
            x1r = a[j] - a[j1];
            x1i = a[j + 1] - a[j1 + 1];
            x2r = a[j2] + a[j3];
            x2i = a[j2 + 1] + a[j3 + 1];
            x3r = a[j2] - a[j3];
            x3i = a[j2 + 1] - a[j3 + 1];
            a[j] = x0r + x2r;
            a[j + 1] = x0i + x2i;
            x0r -= x2r;
            x0i -= x2i;
            a[j2] = wk2r * x0r - wk2i * x0i;
            a[j2 + 1] = wk2r * x0i + wk2i * x0r;
            x0r = x1r - x3i;
            x0i = x1i + x3r;
            a[j1] = wk1r * x0r - wk1i * x0i;
            a[j1 + 1] = wk1r * x0i + wk1i * x0r;
            x0r = x1r + x3i;
            x0i = x1i - x3r;
            a[j3] = wk3r * x0r - wk3i * x0i;
            a[j3 + 1] = wk3r * x0i + wk3i * x0r;
        }
        wk1r = w[k2 + 2];
        wk1i = w[k2 + 3];
        wk3r = wk1r - 2 * wk2r * wk1i;
        wk3i = 2 * wk2r * wk1r - wk1i;
        for (j = k + m; j < l + (k + m); j += 2) {
            j1 = j + l;
            j2 = j1 + l;
            j3 = j2 + l;
            x0r = a[j] + a[j1];
            x0i = a[j + 1] + a[j1 + 1];
            x1r = a[j] - a[j1];
            x1i = a[j + 1] - a[j1 + 1];
            x2r = a[j2] + a[j3];
            x2i = a[j2 + 1] + a[j3 + 1];
            x3r = a[j2] - a[j3];
            x3i = a[j2 + 1] - a[j3 + 1];
            a[j] = x0r + x2r;
            a[j + 1] = x0i + x2i;
            x0r -= x2r;
            x0i -= x2i;
            a[j2] = -wk2i * x0r - wk2r * x0i;
            a[j2 + 1] = -wk2i * x0i + wk2r * x0r;
            x0r = x1r - x3i;
            x0i = x1i + x3r;
            a[j1] = wk1r * x0r - wk1i * x0i;
            a[j1 + 1] = wk1r * x0i + wk1i * x0r;
            x0r = x1r + x3i;
            x0i = x1i - x3r;
            a[j3] = wk3r * x0r - wk3i * x0i;
            a[j3 + 1] = wk3r * x0i + wk3i * x0r;
        }
    }
}

void rftfsub(long n, double *a, int nc, double *c) {
    long j, k, kk, ks, m;
    double wkr, wki, xr, xi, yr, yi;

    m = n >> 1;
    ks = 2 * nc / m;
    kk = 0;
    for (j = 2; j < m; j += 2) {
        k = n - j;
        kk += ks;
        wkr = 0.5 - c[nc - kk];
        wki = c[kk];
        xr = a[j] - a[k];
        xi = a[j + 1] + a[k + 1];
        yr = wkr * xr - wki * xi;
        yi = wkr * xi + wki * xr;
        a[j] -= yr;
        a[j + 1] -= yi;
        a[k] += yr;
        a[k + 1] -= yi;
    }
}

void rftbsub(long n, double *a, int nc, double *c) {
    long j, k, kk, ks, m;
    double wkr, wki, xr, xi, yr, yi;

    a[1] = -a[1];
    m = n >> 1;
    ks = 2 * nc / m;
    kk = 0;
    for (j = 2; j < m; j += 2) {
        k = n - j;
        kk += ks;
        wkr = 0.5 - c[nc - kk];
        wki = c[kk];
        xr = a[j] - a[k];
        xi = a[j + 1] + a[k + 1];
        yr = wkr * xr + wki * xi;
        yi = wkr * xi - wki * xr;
        a[j] -= yr;
        a[j + 1] = yi - a[j + 1];
        a[k] += yr;
        a[k + 1] = yi - a[k + 1];
    }
    a[m + 1] = -a[m + 1];
}

void dctsub(long n, double *a, int nc, double *c) {
    long j, k, kk, ks, m;
    double wkr, wki, xr;

    m = n >> 1;
    ks = nc / n;
    kk = 0;
    for (j = 1; j < m; j++) {
        k = n - j;
        kk += ks;
        wkr = c[kk] - c[nc - kk];
        wki = c[kk] + c[nc - kk];
        xr = wki * a[j] - wkr * a[k];
        a[j] = wkr * a[j] + wki * a[k];
        a[k] = xr;
    }
    a[m] *= c[0];
}

void dstsub(long n, double *a, int nc, double *c) {
    long j, k, kk, ks, m;
    double wkr, wki, xr;

    m = n >> 1;
    ks = nc / n;
    kk = 0;
    for (j = 1; j < m; j++) {
        k = n - j;
        kk += ks;
        wkr = c[kk] - c[nc - kk];
        wki = c[kk] + c[nc - kk];
        xr = wki * a[k] - wkr * a[j];
        a[k] = wkr * a[k] + wki * a[j];
        a[j] = xr;
    }
    a[m] *= c[0];
}

/* sacio.c  2006.4 FUJI Nobuaki
   NPTS will be a power of 2. in readsac
*/

double *readsac(char *sacfile, struct sacheader *hv, long nptsmodifier) {
    FILE *fp;
    double *sacdata;
    fp = fopen(sacfile, "r");
    float ftmp[70];

    char tmpchar[9], tmpchar2[17];
    char ch;
    long i;
    float tmpfloat;

    fread(ftmp, 4, 70, fp);
    hv->DELTA = (double)ftmp[0];
    hv->DEPMIN = (double)ftmp[1];
    hv->DEPMAX = (double)ftmp[2];
    hv->SCALE = (double)ftmp[3];
    hv->ODELTA = (double)ftmp[4];
    hv->B = (double)ftmp[5];
    hv->E = (double)ftmp[6];
    hv->O = (double)ftmp[7];
    hv->A = (double)ftmp[8];
    hv->H9 = (double)ftmp[9];
    hv->T0 = (double)ftmp[10];
    hv->T1 = (double)ftmp[11];
    hv->T2 = (double)ftmp[12];
    hv->T3 = (double)ftmp[13];
    hv->T4 = (double)ftmp[14];
    hv->T5 = (double)ftmp[15];
    hv->T6 = (double)ftmp[16];
    hv->T7 = (double)ftmp[17];
    hv->T8 = (double)ftmp[18];
    hv->T9 = (double)ftmp[19];
    hv->F = (double)ftmp[20];
    hv->RESP0 = (double)ftmp[21];
    hv->RESP1 = (double)ftmp[22];
    hv->RESP2 = (double)ftmp[23];
    hv->RESP3 = (double)ftmp[24];
    hv->RESP4 = (double)ftmp[25];
    hv->RESP5 = (double)ftmp[26];
    hv->RESP6 = (double)ftmp[27];
    hv->RESP7 = (double)ftmp[28];
    hv->RESP8 = (double)ftmp[29];
    hv->RESP9 = (double)ftmp[30];
    hv->STLA = (double)ftmp[31];
    hv->STLO = (double)ftmp[32];
    hv->STEL = (double)ftmp[33];
    hv->STDP = (double)ftmp[34];
    hv->EVLA = (double)ftmp[35];
    hv->EVLO = (double)ftmp[36];
    hv->EVEL = (double)ftmp[37];
    hv->EVDP = (double)ftmp[38];
    hv->H39 = (double)ftmp[39];
    hv->USER0 = (double)ftmp[40];
    hv->USER1 = (double)ftmp[41];
    hv->USER2 = (double)ftmp[42];
    hv->USER3 = (double)ftmp[43];
    hv->USER4 = (double)ftmp[44];
    hv->USER5 = (double)ftmp[45];
    hv->USER6 = (double)ftmp[46];
    hv->USER7 = (double)ftmp[47];
    hv->USER8 = (double)ftmp[48];
    hv->USER9 = (double)ftmp[49];
    hv->DIST = (double)ftmp[50];
    hv->AZ = (double)ftmp[51];
    hv->BAZ = (double)ftmp[52];
    hv->GCARC = (double)ftmp[53];
    hv->H54 = (double)ftmp[54];
    hv->H55 = (double)ftmp[55];
    hv->DEPMEN = (double)ftmp[56];
    hv->CMPAZ = (double)ftmp[57];
    hv->CMPINC = (double)ftmp[58];
    hv->H59 = (double)ftmp[59];
    hv->H60 = (double)ftmp[60];
    hv->H61 = (double)ftmp[61];
    hv->H62 = (double)ftmp[62];
    hv->H63 = (double)ftmp[63];
    hv->H64 = (double)ftmp[64];
    hv->H65 = (double)ftmp[65];
    hv->H66 = (double)ftmp[66];
    hv->H67 = (double)ftmp[67];
    hv->H68 = (double)ftmp[68];
    hv->H69 = (double)ftmp[69];
    fread(&hv->NZYEAR, 4, 1, fp);
    fread(&hv->NZJDAY, 4, 1, fp);
    fread(&hv->NZHOUR, 4, 1, fp);
    fread(&hv->NZMIN, 4, 1, fp);
    fread(&hv->NZSEC, 4, 1, fp);
    fread(&hv->NZMSEC, 4, 1, fp);
    fread(&hv->NVHDR, 4, 1, fp);
    fread(&hv->H77, 4, 1, fp);
    fread(&hv->H78, 4, 1, fp);
    fread(&hv->NPTS, 4, 1, fp);
    fread(&hv->H80, 4, 1, fp);
    fread(&hv->H81, 4, 1, fp);
    fread(&hv->H82, 4, 1, fp);
    fread(&hv->H83, 4, 1, fp);
    fread(&hv->H84, 4, 1, fp);
    fread(&hv->IFTYPE, 4, 1, fp);
    fread(&hv->IDEP, 4, 1, fp);
    fread(&hv->IZTYPE, 4, 1, fp);
    fread(&hv->H88, 4, 1, fp);
    fread(&hv->H89, 4, 1, fp);
    fread(&hv->ISTREG, 4, 1, fp);
    fread(&hv->IEVREG, 4, 1, fp);
    fread(&hv->IEVTYP, 4, 1, fp);
    fread(&hv->IQUAL, 4, 1, fp);
    fread(&hv->ISYNTH, 4, 1, fp);
    fread(&hv->H95, 4, 1, fp);
    fread(&hv->H96, 4, 1, fp);
    fread(&hv->H97, 4, 1, fp);
    fread(&hv->H98, 4, 1, fp);
    fread(&hv->H99, 4, 1, fp);
    fread(&hv->H100, 4, 1, fp);
    fread(&hv->H101, 4, 1, fp);
    fread(&hv->H102, 4, 1, fp);
    fread(&hv->H103, 4, 1, fp);
    fread(&hv->H104, 4, 1, fp);
    fread(&hv->LEVEN, 4, 1, fp);
    fread(&hv->LPSPOL, 4, 1, fp);
    fread(&hv->LOVROK, 4, 1, fp);
    fread(&hv->LCALDA, 4, 1, fp);
    fread(&hv->H109, 4, 1, fp);

    for (i = 0; i < 8; i++) {
        fread(&ch, 1, 1, fp);
        tmpchar[i] = ch;
    }
    tmpchar[8] = '\0';
    strcpy(hv->KSTNM, tmpchar);

    for (i = 0; i < 16; i++) {
        fread(&ch, 1, 1, fp);
        tmpchar2[i] = ch;
    }
    tmpchar2[16] = '\0';
    strcpy(hv->KEVNM, tmpchar2);

    for (i = 0; i < 8; i++) {
        fread(&ch, 1, 1, fp);
        tmpchar[i] = ch;
    }
    tmpchar[8] = '\0';
    strcpy(hv->KHOLE, tmpchar);

    for (i = 0; i < 8; i++) {
        fread(&ch, 1, 1, fp);
        tmpchar[i] = ch;
    }
    tmpchar[8] = '\0';
    strcpy(hv->KO, tmpchar);

    for (i = 0; i < 8; i++) {
        fread(&ch, 1, 1, fp);
        tmpchar[i] = ch;
    }
    tmpchar[8] = '\0';
    strcpy(hv->KA, tmpchar);

    for (i = 0; i < 8; i++) {
        fread(&ch, 1, 1, fp);
        tmpchar[i] = ch;
    }
    tmpchar[8] = '\0';
    strcpy(hv->KT0, tmpchar);

    for (i = 0; i < 8; i++) {
        fread(&ch, 1, 1, fp);
        tmpchar[i] = ch;
    }
    tmpchar[8] = '\0';
    strcpy(hv->KT1, tmpchar);

    for (i = 0; i < 8; i++) {
        fread(&ch, 1, 1, fp);
        tmpchar[i] = ch;
    }
    tmpchar[8] = '\0';
    strcpy(hv->KT2, tmpchar);

    for (i = 0; i < 8; i++) {
        fread(&ch, 1, 1, fp);
        tmpchar[i] = ch;
    }
    tmpchar[8] = '\0';
    strcpy(hv->KT3, tmpchar);

    for (i = 0; i < 8; i++) {
        fread(&ch, 1, 1, fp);
        tmpchar[i] = ch;
    }
    tmpchar[8] = '\0';
    strcpy(hv->KT4, tmpchar);

    for (i = 0; i < 8; i++) {
        fread(&ch, 1, 1, fp);
        tmpchar[i] = ch;
    }
    tmpchar[8] = '\0';
    strcpy(hv->KT5, tmpchar);

    for (i = 0; i < 8; i++) {
        fread(&ch, 1, 1, fp);
        tmpchar[i] = ch;
    }
    tmpchar[8] = '\0';
    strcpy(hv->KT6, tmpchar);

    for (i = 0; i < 8; i++) {
        fread(&ch, 1, 1, fp);
        tmpchar[i] = ch;
    }
    tmpchar[8] = '\0';
    strcpy(hv->KT7, tmpchar);

    for (i = 0; i < 8; i++) {
        fread(&ch, 1, 1, fp);
        tmpchar[i] = ch;
    }
    tmpchar[8] = '\0';
    strcpy(hv->KT8, tmpchar);

    for (i = 0; i < 8; i++) {
        fread(&ch, 1, 1, fp);
        tmpchar[i] = ch;
    }
    tmpchar[8] = '\0';
    strcpy(hv->KT9, tmpchar);

    for (i = 0; i < 8; i++) {
        fread(&ch, 1, 1, fp);
        tmpchar[i] = ch;
    }
    tmpchar[8] = '\0';
    strcpy(hv->KF, tmpchar);

    for (i = 0; i < 8; i++) {
        fread(&ch, 1, 1, fp);
        tmpchar[i] = ch;
    }
    tmpchar[8] = '\0';
    strcpy(hv->KUSER0, tmpchar);

    for (i = 0; i < 8; i++) {
        fread(&ch, 1, 1, fp);
        tmpchar[i] = ch;
    }
    tmpchar[8] = '\0';
    strcpy(hv->KUSER1, tmpchar);

    for (i = 0; i < 8; i++) {
        fread(&ch, 1, 1, fp);
        tmpchar[i] = ch;
    }
    tmpchar[8] = '\0';
    strcpy(hv->KUSER2, tmpchar);

    for (i = 0; i < 8; i++) {
        fread(&ch, 1, 1, fp);
        tmpchar[i] = ch;
    }
    tmpchar[8] = '\0';
    strcpy(hv->KCMPNM, tmpchar);

    for (i = 0; i < 8; i++) {
        fread(&ch, 1, 1, fp);
        tmpchar[i] = ch;
    }
    tmpchar[8] = '\0';
    strcpy(hv->KNETWK, tmpchar);

    for (i = 0; i < 8; i++) {
        fread(&ch, 1, 1, fp);
        tmpchar[i] = ch;
    }
    tmpchar[8] = '\0';
    strcpy(hv->KDATRD, tmpchar);

    for (i = 0; i < 8; i++) {
        fread(&ch, 1, 1, fp);
        tmpchar[i] = ch;
    }
    tmpchar[8] = '\0';
    strcpy(hv->KINST, tmpchar);

    if (nptsmodifier == 1) {
        while (nptsmodifier < hv->NPTS) nptsmodifier *= 2;

        sacdata = calloc((size_t)hv->NPTS, sizeof(double));

        if (sacdata == NULL) {
            printf("We could not obtain enough memory for sacdata.\n");
            exit(1);
        }

        for (i = 0; i < hv->NPTS; i++) {
            fread(&tmpfloat, 4, 1, fp);
            sacdata[i] = (double)tmpfloat;
        }

        if (hv->NPTS < nptsmodifier) {
            for (i = hv->NPTS; i < nptsmodifier; i++) {
                sacdata[i] = 0.0;
            }
            hv->NPTS = i;
        }
    } else {
        sacdata = calloc((size_t)hv->NPTS, sizeof(double));

        if (sacdata == NULL) {
            printf("We could not obtain enough memory for sacdata.\n");
            exit(1);
        }

        for (i = 0; i < hv->NPTS; i++) {
            fread(&tmpfloat, 4, 1, fp);
            sacdata[i] = (double)tmpfloat;
        }
    }
    fclose(fp);

    return (sacdata);
}

int newsacheader(struct sacheader *hv) {
    /* Thank you UCHIDE-san */

    hv->DELTA = -12345.0;
    hv->DEPMIN = -12345.0;
    hv->DEPMAX = -12345.0;
    hv->SCALE = -12345.0;
    hv->ODELTA = -12345.0;
    hv->B = -12345.0;
    hv->E = -12345.0;
    hv->O = -12345.0;
    hv->A = -12345.0;
    hv->H9 = -12345.0;
    hv->T0 = -12345.0;
    hv->T1 = -12345.0;
    hv->T2 = -12345.0;
    hv->T3 = -12345.0;
    hv->T4 = -12345.0;
    hv->T5 = -12345.0;
    hv->T6 = -12345.0;
    hv->T7 = -12345.0;
    hv->T8 = -12345.0;
    hv->T9 = -12345.0;
    hv->F = -12345.0;
    hv->RESP0 = -12345.0;
    hv->RESP1 = -12345.0;
    hv->RESP2 = -12345.0;
    hv->RESP3 = -12345.0;
    hv->RESP4 = -12345.0;
    hv->RESP5 = -12345.0;
    hv->RESP6 = -12345.0;
    hv->RESP7 = -12345.0;
    hv->RESP8 = -12345.0;
    hv->RESP9 = -12345.0;
    hv->STLA = -12345.0;
    hv->STLO = -12345.0;
    hv->STEL = -12345.0;
    hv->STDP = -12345.0;
    hv->EVLA = -12345.0;
    hv->EVLO = -12345.0;
    hv->EVEL = -12345.0;
    hv->EVDP = -12345.0;
    hv->H39 = -12345.0;
    hv->USER0 = -12345.0;
    hv->USER1 = -12345.0;
    hv->USER2 = -12345.0;
    hv->USER3 = -12345.0;
    hv->USER4 = -12345.0;
    hv->USER5 = -12345.0;
    hv->USER6 = -12345.0;
    hv->USER7 = -12345.0;
    hv->USER8 = -12345.0;
    hv->USER9 = -12345.0;
    hv->DIST = -12345.0;
    hv->AZ = -12345.0;
    hv->BAZ = -12345.0;
    hv->GCARC = -12345.0;
    ;
    hv->H54 = -12345.0;
    hv->H55 = -12345.0;
    hv->DEPMEN = -12345.0;
    hv->CMPAZ = -12345.0;
    hv->CMPINC = -12345.0;
    hv->H59 = -12345.0;
    hv->H60 = -12345.0;
    hv->H61 = -12345.0;
    hv->H62 = -12345.0;
    hv->H63 = -12345.0;
    hv->H64 = -12345.0;
    hv->H65 = -12345.0;
    hv->H66 = -12345.0;
    hv->H67 = -12345.0;
    hv->H68 = -12345.0;
    hv->H69 = -12345.0;
    hv->NZYEAR = -12345;
    hv->NZJDAY = -12345;
    hv->NZHOUR = -12345;
    hv->NZMIN = -12345;
    hv->NZSEC = -12345;
    hv->NZMSEC = -12345;
    hv->NVHDR = 6;
    hv->H77 = -12345;
    hv->H78 = -12345;
    hv->NPTS = -12345;
    hv->H80 = -12345;
    hv->H81 = -12345;
    hv->H82 = -12345;
    hv->H83 = -12345;
    hv->H84 = -12345;
    hv->IFTYPE = 1;
    hv->IDEP = 7;
    hv->IZTYPE = -12345;
    hv->H88 = -12345;
    hv->H89 = -12345;
    hv->ISTREG = -12345;
    hv->IEVREG = -12345;
    hv->IEVTYP = -12345;
    hv->IQUAL = -12345;
    hv->ISYNTH = -12345;
    hv->H95 = -12345;
    hv->H96 = -12345;
    hv->H97 = -12345;
    hv->H98 = -12345;
    hv->H99 = -12345;
    hv->H100 = -12345;
    hv->H101 = -12345;
    hv->H102 = -12345;
    hv->H103 = -12345;
    hv->H104 = -12345;
    hv->LEVEN = 1;
    hv->LPSPOL = 0;
    hv->LOVROK = 1;
    hv->LCALDA = 1;
    hv->H109 = -12345;
    strcpy(hv->KSTNM, "-12345");
    strcpy(hv->KEVNM, "-12345");
    strcpy(hv->KHOLE, "-12345");
    strcpy(hv->KO, "-12345");
    strcpy(hv->KA, "-12345");
    strcpy(hv->KT0, "-12345");
    strcpy(hv->KT1, "-12345");
    strcpy(hv->KT2, "-12345");
    strcpy(hv->KT3, "-12345");
    strcpy(hv->KT4, "-12345");
    strcpy(hv->KT5, "-12345");
    strcpy(hv->KT6, "-12345");
    strcpy(hv->KT7, "-12345");
    strcpy(hv->KT8, "-12345");
    strcpy(hv->KT9, "-12345");
    strcpy(hv->KF, "-12345");
    strcpy(hv->KUSER0, "-12345");
    strcpy(hv->KUSER1, "-12345");
    strcpy(hv->KUSER2, "-12345");
    strcpy(hv->KCMPNM, "-12345");
    strcpy(hv->KNETWK, "-12345");
    strcpy(hv->KDATRD, "-12345");
    strcpy(hv->KINST, "-12345");

    return 0;
}

int writeasciisac(char *sacfile, struct sacheader *hv, double *sacdata) {
    FILE *fp;
    fp = (fopen(sacfile, "w"));
    long i;
    for (i = 0; i < hv->NPTS; i++) {
        fprintf(fp, "%e %e\n", hv->DELTA * (double)i, sacdata[i]);
    }
    fclose(fp);
    return 0;
}

int writesac(char *sacfile, struct sacheader *hv, double *sacdata) {
    FILE *fp;
    fp = fopen(sacfile, "w");
    float ftmp[70], tmpfloat;
    long i;
    char ch;
    char tmpchar[9], tmpchar2[17];

    ftmp[0] = (float)hv->DELTA;
    ftmp[1] = (float)hv->DEPMIN;
    ftmp[2] = (float)hv->DEPMAX;
    ftmp[3] = (float)hv->SCALE;
    ftmp[4] = (float)hv->ODELTA;
    ftmp[5] = (float)hv->B;
    ftmp[6] = (float)hv->E;
    ftmp[7] = (float)hv->O;
    ftmp[8] = (float)hv->A;
    ftmp[9] = (float)hv->H9;
    ftmp[10] = (float)hv->T0;
    ftmp[11] = (float)hv->T1;
    ftmp[12] = (float)hv->T2;
    ftmp[13] = (float)hv->T3;
    ftmp[14] = (float)hv->T4;
    ftmp[15] = (float)hv->T5;
    ftmp[16] = (float)hv->T6;
    ftmp[17] = (float)hv->T7;
    ftmp[18] = (float)hv->T8;
    ftmp[19] = (float)hv->T9;
    ftmp[20] = (float)hv->F;
    ftmp[21] = (float)hv->RESP0;
    ftmp[22] = (float)hv->RESP1;
    ftmp[23] = (float)hv->RESP2;
    ftmp[24] = (float)hv->RESP3;
    ftmp[25] = (float)hv->RESP4;
    ftmp[26] = (float)hv->RESP5;
    ftmp[27] = (float)hv->RESP6;
    ftmp[28] = (float)hv->RESP7;
    ftmp[29] = (float)hv->RESP8;
    ftmp[30] = (float)hv->RESP9;
    ftmp[31] = (float)hv->STLA;
    ftmp[32] = (float)hv->STLO;
    ftmp[33] = (float)hv->STEL;
    ftmp[34] = (float)hv->STDP;
    ftmp[35] = (float)hv->EVLA;
    ftmp[36] = (float)hv->EVLO;
    ftmp[37] = (float)hv->EVEL;
    ftmp[38] = (float)hv->EVDP;
    ftmp[39] = (float)hv->H39;
    ftmp[40] = (float)hv->USER0;
    ftmp[41] = (float)hv->USER1;
    ftmp[42] = (float)hv->USER2;
    ftmp[43] = (float)hv->USER3;
    ftmp[44] = (float)hv->USER4;
    ftmp[45] = (float)hv->USER5;
    ftmp[46] = (float)hv->USER6;
    ftmp[47] = (float)hv->USER7;
    ftmp[48] = (float)hv->USER8;
    ftmp[49] = (float)hv->USER9;
    ftmp[50] = (float)hv->DIST;
    ftmp[51] = (float)hv->AZ;
    ftmp[52] = (float)hv->BAZ;
    ftmp[53] = (float)hv->GCARC;
    ftmp[54] = (float)hv->H54;
    ftmp[55] = (float)hv->H55;
    ftmp[56] = (float)hv->DEPMEN;
    ftmp[57] = (float)hv->CMPAZ;
    ftmp[58] = (float)hv->CMPINC;
    ftmp[59] = (float)hv->H59;
    ftmp[60] = (float)hv->H60;
    ftmp[61] = (float)hv->H61;
    ftmp[62] = (float)hv->H62;
    ftmp[63] = (float)hv->H63;
    ftmp[64] = (float)hv->H64;
    ftmp[65] = (float)hv->H65;
    ftmp[66] = (float)hv->H66;
    ftmp[67] = (float)hv->H67;
    ftmp[68] = (float)hv->H68;
    ftmp[69] = (float)hv->H69;
    fwrite(ftmp, 4, 70, fp);
    fwrite(&hv->NZYEAR, 4, 1, fp);
    fwrite(&hv->NZJDAY, 4, 1, fp);
    fwrite(&hv->NZHOUR, 4, 1, fp);
    fwrite(&hv->NZMIN, 4, 1, fp);
    fwrite(&hv->NZSEC, 4, 1, fp);
    fwrite(&hv->NZMSEC, 4, 1, fp);
    fwrite(&hv->NVHDR, 4, 1, fp);
    fwrite(&hv->H77, 4, 1, fp);
    fwrite(&hv->H78, 4, 1, fp);
    fwrite(&hv->NPTS, 4, 1, fp);
    fwrite(&hv->H80, 4, 1, fp);
    fwrite(&hv->H81, 4, 1, fp);
    fwrite(&hv->H82, 4, 1, fp);
    fwrite(&hv->H83, 4, 1, fp);
    fwrite(&hv->H84, 4, 1, fp);
    fwrite(&hv->IFTYPE, 4, 1, fp);
    fwrite(&hv->IDEP, 4, 1, fp);
    fwrite(&hv->IZTYPE, 4, 1, fp);
    fwrite(&hv->H88, 4, 1, fp);
    fwrite(&hv->H89, 4, 1, fp);
    fwrite(&hv->ISTREG, 4, 1, fp);
    fwrite(&hv->IEVREG, 4, 1, fp);
    fwrite(&hv->IEVTYP, 4, 1, fp);
    fwrite(&hv->IQUAL, 4, 1, fp);
    fwrite(&hv->ISYNTH, 4, 1, fp);
    fwrite(&hv->H95, 4, 1, fp);
    fwrite(&hv->H96, 4, 1, fp);
    fwrite(&hv->H97, 4, 1, fp);
    fwrite(&hv->H98, 4, 1, fp);
    fwrite(&hv->H99, 4, 1, fp);
    fwrite(&hv->H100, 4, 1, fp);
    fwrite(&hv->H101, 4, 1, fp);
    fwrite(&hv->H102, 4, 1, fp);
    fwrite(&hv->H103, 4, 1, fp);
    fwrite(&hv->H104, 4, 1, fp);
    fwrite(&hv->LEVEN, 4, 1, fp);
    fwrite(&hv->LPSPOL, 4, 1, fp);
    fwrite(&hv->LOVROK, 4, 1, fp);
    fwrite(&hv->LCALDA, 4, 1, fp);
    fwrite(&hv->H109, 4, 1, fp);

    strcpy(tmpchar, hv->KSTNM);
    for (i = 0; i < 8; i++) {
        ch = tmpchar[i];
        fwrite(&ch, 1, 1, fp);
    }

    strcpy(tmpchar2, hv->KEVNM);
    for (i = 0; i < 16; i++) {
        ch = tmpchar2[i];
        fwrite(&ch, 1, 1, fp);
    }

    strcpy(tmpchar, hv->KHOLE);
    for (i = 0; i < 8; i++) {
        ch = tmpchar[i];
        fwrite(&ch, 1, 1, fp);
    }

    strcpy(tmpchar, hv->KO);
    for (i = 0; i < 8; i++) {
        ch = tmpchar[i];
        fwrite(&ch, 1, 1, fp);
    }

    strcpy(tmpchar, hv->KA);
    for (i = 0; i < 8; i++) {
        ch = tmpchar[i];
        fwrite(&ch, 1, 1, fp);
    }

    strcpy(tmpchar, hv->KT0);
    for (i = 0; i < 8; i++) {
        ch = tmpchar[i];
        fwrite(&ch, 1, 1, fp);
    }

    strcpy(tmpchar, hv->KT1);
    for (i = 0; i < 8; i++) {
        ch = tmpchar[i];
        fwrite(&ch, 1, 1, fp);
        tmpchar[i] = ch;
    }

    strcpy(tmpchar, hv->KT2);
    for (i = 0; i < 8; i++) {
        ch = tmpchar[i];
        fwrite(&ch, 1, 1, fp);
    }

    strcpy(tmpchar, hv->KT3);
    for (i = 0; i < 8; i++) {
        ch = tmpchar[i];
        fwrite(&ch, 1, 1, fp);
    }

    strcpy(tmpchar, hv->KT4);
    for (i = 0; i < 8; i++) {
        ch = tmpchar[i];
        fwrite(&ch, 1, 1, fp);
    }

    strcpy(tmpchar, hv->KT5);
    for (i = 0; i < 8; i++) {
        ch = tmpchar[i];
        fwrite(&ch, 1, 1, fp);
    }

    strcpy(tmpchar, hv->KT6);
    for (i = 0; i < 8; i++) {
        ch = tmpchar[i];
        fwrite(&ch, 1, 1, fp);
    }

    strcpy(tmpchar, hv->KT7);
    for (i = 0; i < 8; i++) {
        ch = tmpchar[i];
        fwrite(&ch, 1, 1, fp);
    }

    strcpy(tmpchar, hv->KT8);
    for (i = 0; i < 8; i++) {
        ch = tmpchar[i];
        fwrite(&ch, 1, 1, fp);
    }

    strcpy(tmpchar, hv->KT9);
    for (i = 0; i < 8; i++) {
        ch = tmpchar[i];
        fwrite(&ch, 1, 1, fp);
    }

    strcpy(tmpchar, hv->KF);
    for (i = 0; i < 8; i++) {
        ch = tmpchar[i];
        fwrite(&ch, 1, 1, fp);
    }

    strcpy(tmpchar, hv->KUSER0);
    for (i = 0; i < 8; i++) {
        ch = tmpchar[i];
        fwrite(&ch, 1, 1, fp);
    }

    strcpy(tmpchar, hv->KUSER1);
    for (i = 0; i < 8; i++) {
        ch = tmpchar[i];
        fwrite(&ch, 1, 1, fp);
    }

    strcpy(tmpchar, hv->KUSER2);
    for (i = 0; i < 8; i++) {
        ch = tmpchar[i];
        fwrite(&ch, 1, 1, fp);
    }

    strcpy(tmpchar, hv->KCMPNM);
    for (i = 0; i < 8; i++) {
        ch = tmpchar[i];
        fwrite(&ch, 1, 1, fp);
        tmpchar[i] = ch;
    }

    strcpy(tmpchar, hv->KNETWK);
    for (i = 0; i < 8; i++) {
        ch = tmpchar[i];
        fwrite(&ch, 1, 1, fp);
    }

    strcpy(tmpchar, hv->KDATRD);
    for (i = 0; i < 8; i++) {
        ch = tmpchar[i];
        fwrite(&ch, 1, 1, fp);
    }

    strcpy(tmpchar, hv->KINST);
    for (i = 0; i < 8; i++) {
        ch = tmpchar[i];
        fwrite(&ch, 1, 1, fp);
    }

    for (i = 0; i < hv->NPTS; i++) {
        tmpfloat = (float)sacdata[i];
        fwrite(&tmpfloat, 4, 1, fp);
    }

    return (0);
}
