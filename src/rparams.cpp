
char filename[256];
char radialfilename[256];
unsigned short key[3];
int num_steps = 1000;
int num_divisions = 2;
float bin_width = 1.0;

void printhelp() {
  printf("\nformat:  radial <sequence filename> <radial histogram> [extra arguments]\n");
  printf("\nExtra arguments [defaults]:\n\n");
  printf("-n integer[1000]  Changes the number of samples to take\n");
  printf("-b number[1.0]    Sets the bin width for radial density distribution\n");
  printf("-k number[0]		Sets the value of the random number key used\n");
  printf("Example: ./radial seq400 radial.dat -n 10000 -k 1 -b 5.0\n");
}

void parsecommandline(int argc, char **argv) {

  if (argc == 1) {
    printhelp();
    exit(0);
  }

  strcpy(filename, argv[1]);
  strcpy(radialfilename, argv[2]);

  for (int i = 3; i < argc; i++) {

    if (strcmp(argv[i], "-n") == 0) {
      sscanf(argv[i+1], "%d", &num_steps);
      i++;
    } else
    if (strcmp(argv[i], "-k") == 0) {
      sscanf(argv[i+1], "%d", &key[0]);
      key[1] = 0;
      key[2] = 0;
      printf("Using key %d 0 0\n", key[0]);
      i++;
    } else
    if (strcmp(argv[i], "-b") == 0) {
	  sscanf(argv[i+1], "%f", &bin_width);
      i++;
    } else {
      printf("Garbage encountered on command line: %s\n\n", argv[i]);
      exit(0);
    }	       
  }
  
}
