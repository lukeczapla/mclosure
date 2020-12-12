
char filename[256];
FLOAT bound[20][6];
FLOAT epsilon = 10.0;
FLOAT beta = 1.0;
FLOAT radi = 30;
FLOAT gam = 0.866;
FLOAT twi = 30;
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
int persistence = 0;
int regen = 0;
int reduced = 0;
long int key;
int num_HU = 0;
int make_pictures = 0;
int make_spictures = 0;

void printhelp() {
  printf("format:  closuremc <sequence filename> [extra arguments]\n");
  printf("Extra arguments [defaults]:\n");
  printf("-e number[10.0]     Changes the epsilon (cube length of \"bin\")\n");
  printf("-P                  Calculate persistence length\n");
  printf("-n integer[1000]    Changes the number of samples to take\n");
  printf("-b 6 numbers        Adds an end-to-end theta/rho boundary condition\n");
  printf("-s integer[0]       Number of confs to save to memory for closure-disk\n");
  printf("-B number[1.0]      Sets the value of beta\n");
  printf("-R                  Regenerate the system with the pair file\n");
  printf("-pic                Generates pictures (used with -R)\n");
  printf("-rpic               Generates supercoil pictures (used with -R)\n");
#ifdef _HU
  printf("-HU number[0.0129]  HU binding probability (greater than 0)\n");
  printf("-HUpos filename     Creates a closure HU position histogram\n");
  printf("-HUdis filename     Creates a closure HU distance histogram\n");
  printf("-HUn filename       Creates a closure N_HU histogram\n");
#endif
#ifdef _HU2
  printf("-H number[0]        Number of HUs to place per segment\n");
  printf("-HUpos filename     Creates a closure HU position histogram\n");
  printf("-HUdis filename     Creates a closure HU distance histogram\n");
  printf("-HUn filename       Creates a closure N_HU histogram\n");
#endif
  printf("-radi number[30]    (Angs) sets radial bound r\n");
  printf("-gam number[0.86]   sets cos(gamma) > bound \n");
  printf("-twi number[30]     (deg) sets abs(twist) bound \n");
  printf("-k number[0]        sets the srand48 key\n");
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
    if (strcmp(argv[i], "-n") == 0) {
      sscanf(argv[i+1], "%d", &num_steps);
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
    if (strcmp(argv[i], "-d") == 0) {
      sscanf(argv[i+1], "%d", &num_divisions);
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
    if (strcmp(argv[i], "-s") == 0) {
      sscanf(argv[i+1], "%d", &num_saved);
      i++;
    } else 
    if (strcmp(argv[i], "-k") == 0) {
      sscanf(argv[i+1], "%ld", &key);
      printf("Using key %ld\n", key);
      i++;
    } else {
      printf("Garbage encountered on command line: %s\n\n", argv[i]);
      exit(0);
    }	       
  }
  
}
