
char filename[256];
FLOAT bound[20][6];
FLOAT epsilon = 10.0;
FLOAT beta = 1.0;
FLOAT radi = 30;
FLOAT gam = 0.866;
FLOAT twi = 30;
FLOAT sc_width = 0.9;
FLOAT bin_width = 0.05;
#ifdef _HU
FLOAT PHU = 0.0129;
#else
FLOAT PHU = 0.0;
#endif
char positionH[30];
char spacingH[30];
char nhuH[30];
int num_steps = 10000;
int num_divisions = 2;
int num_saved = 0;
int j = 0;
int restrictends = 1;
int persistence = 0;
int regen = 0;
int reduced = 0;
int overlap = 0;
unsigned short key[3];
int osamples = 0;
int neglect = 0;
int num_HU = 0;
int make_pictures = 0;
int make_spictures = 0;
int himem = 0;

void printhelp() {
  printf("\nformat:  closuremc <sequence filename> [extra arguments]\n");
  printf("\nExtra arguments [defaults]:\n\n");
  printf("-P                Calculate persistence length and force extension data\n");
  printf("-r                Reduce the complexity of quadratic combinations\n");
  printf("-n integer[1000]  Changes the number of samples to take\n");
  printf("-B number[1.0]    Sets the value of beta\n");
  printf("-o number[0]      Sets the number of random combinations to take\n");
  printf("-b 6 numbers      Adds a theta/rho boundary condition (default=[0,0,0,0,0,0])\n\n");
  printf("-pic              Generates pictures of closed molecules\n");
  printf("-rpic             Generates supercoiled pictures, but not normal ones\n");
  printf("-sc number[0.9]   Save default supercoiled molecules with dLk >= number\n");
  printf("-wr number [0.05] Sets the bin size for writhe and twist histograms\n");
#ifdef _HU
  printf("-HU number[0.013] HU binding probability (greater than 0)\n");
  printf("-HUpos filename   Creates a closure HU position histogram\n");
  printf("-HUdis filename   Creates a closure HU distance histogram\n");
  printf("-HUn filename     Creates a closure N_HU histogram\n");
#endif
#ifdef _HU2
  printf("-H number[0]      Number of HUs to place per segment\n");
  printf("-HUpos filename   Creates a closure HU position histogram\n");
  printf("-HUdis filename   Creates a closure HU distance histogram\n");
  printf("-HUn filename     Creates a closure N_HU histogram\n");
#endif
  printf("-overlap          Turn on overlap checking for half segments\n");
  printf("-e number[10.0]   Changes the epsilon (cube length of \"bin\")\n");
  printf("-radi number[30]  Sets radial bounds in Angstroms\n");
  printf("-gam number[0.86] Sets cos(gamma) terminal normal bounds, cos(gam) > N\n");
  printf("-twi number[30]   Sets abs(twist) bound in degrees \n");
  printf("-k number[0]      Sets the erand48_0 key\n");
  printf("-loose            Removes hemispheric restriction (z < 0)\n");
  printf("-M                Use extra memory for smaller calculations (N<5,000,000)\n");
  printf("-no               Neglect overlap altogether\n\n");
  printf("Example: ./closuremc seq105 -r -n 100000 -k 42 -e 5.0 -radi 20.0 -gam 0.98 -twi 11.5\n");
}

void parsecommandline(int argc, char **argv) {

  if (argc == 1) {
    printhelp();
    exit(0);
  }

  strcpy(filename, argv[1]);

  strcpy(positionH, "");
  strcpy(nhuH, "");
  strcpy(spacingH, "");

  for (int i = 2; i < argc; i++) {

    if (strcmp(argv[i], "-e") == 0) {
      sscanf(argv[i+1], FLOAT_ESTR, &epsilon);
      i++;
    } else
	if (strcmp(argv[i], "-overlap") == 0) {
	  overlap = 1;
	} else
    if (strcmp(argv[i], "-loose") == 0) {
	  restrictends = 0;
    } else
	if (strcmp(argv[i], "-o") == 0) {
	   sscanf(argv[i+1], "%d", &osamples);
	   i++;
	} else
    if (strcmp(argv[i], "-M") == 0) {
       himem = 1;
       printf("Using extra memory for high performance\n");
    } else
    if (strcmp(argv[i], "-n") == 0) {
      sscanf(argv[i+1], "%d", &num_steps);
      i++;
    } else
	if (strcmp(argv[i], "-sc") == 0) {
	  sscanf(argv[i+1], "%lf", &sc_width);
	  printf("Setting supercoiling dLk > %lf\n", sc_width);
	  i++;
	} else
	if (strcmp(argv[i], "-wr") == 0) {
	  sscanf(argv[i+1], "%lf", &bin_width);
	  i++;
	} else
    if (strcmp(argv[i], "-R") == 0) {
      regen = 1;
    } else
    if (strcmp(argv[i], "-P") == 0) {
      persistence = 1;
    } else
    if (strcmp(argv[i], "-radi") == 0) {
      sscanf(argv[i+1], FLOAT_ESTR, &radi);
      i++;
    } else
    if (strcmp(argv[i], "-H") == 0) {
      sscanf(argv[i+1], "%d", &num_HU);
      i++;
    } else
    if (strcmp(argv[i], "-pic") == 0) {
	  make_pictures = 1;
    } else
	if (strcmp(argv[i], "-rpic") == 0) {
	  make_spictures = 1;
	} else
    if (strcmp(argv[i], "-gam") == 0) {
      sscanf(argv[i+1], FLOAT_ESTR, &gam);
      i++;
    } else
    if (strcmp(argv[i], "-twi") == 0) {
      sscanf(argv[i+1], FLOAT_ESTR, &twi);
      i++;
    } else 
    if (strcmp(argv[i], "-b") == 0) {
      sscanf(argv[i+1], FLOAT_ESTR, &bound[j][0]);
      sscanf(argv[i+2], FLOAT_ESTR, &bound[j][1]);
      sscanf(argv[i+3], FLOAT_ESTR, &bound[j][2]);
      sscanf(argv[i+4], FLOAT_ESTR, &bound[j][3]);
      sscanf(argv[i+5], FLOAT_ESTR, &bound[j][4]);
      sscanf(argv[i+6], FLOAT_ESTR, &bound[j][5]);
      i += 6;
      j++;
    } else
    if (strcmp(argv[i], "-B") == 0) {
      sscanf(argv[i+1], FLOAT_ESTR, &beta);
      printf("beta = %e\n", beta);
      i++;
    } else
    if (strcmp(argv[i], "-HU") == 0) {
      sscanf(argv[i+1], FLOAT_ESTR, &PHU);
      printf("HU probability = %e\n", PHU);
      i++;
    } else
    if (strcmp(argv[i], "-HUn") == 0) {
      sscanf(argv[i+1], "%s", nhuH);
      printf("HU number histogram file = %s\n", nhuH);
      i++;
    } else
    if (strcmp(argv[i], "-HUpos") == 0) {
      sscanf(argv[i+1], "%s", positionH);
      printf("HU position histogram file = %s\n", positionH);
      i++;
    } else
    if (strcmp(argv[i], "-HUdis") == 0) {
      sscanf(argv[i+1], "%s", spacingH);
      printf("HU spacing histogram file = %s\n", spacingH);
      i++;
    } else      
    if (strcmp(argv[i], "-r") == 0) {
      reduced = 1;
    } else
    if (strcmp(argv[i], "-k") == 0) {
      sscanf(argv[i+1], "%d", &key[0]);
      key[1] = 0;
      key[2] = 0;
      printf("Using key %d 0 0\n", key[0]);
      i++;
    } else
    if (strcmp(argv[i], "-no") == 0) {
      neglect = 1;
      i++;
    } else
    if (strcmp(argv[i], "-s") == 0) {
      sscanf(argv[i+1], "%d", &num_saved);
      i++;
    } else {
      printf("Garbage encountered on command line: %s\n\n", argv[i]);
      exit(0);
    }	       
  }
  
}
