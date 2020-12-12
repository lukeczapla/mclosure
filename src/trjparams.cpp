
int ieh = 0;
int iet = 0;
int isp = 0;
int irh = 0;
int irt = 0;
int igh = 0;
int igt = 0;

int start = 0;
int end = 2000000;

char eh[256];
char et[256];
char sp[256];
char rh[256];
char rt[256];
char gh[256];
char gt[256];
char filename[256];
char source[256];


void printhelp() {
  printf("usage:  calctrj <sequence file> <trajectory filename> [extra arguments]\n");
  printf("-i number [0]         - Sets the initial frame to start calculating\n");
  printf("-f number [end]       - Sets the final frame to calculate\n");
  printf("-e <filename>         - Writes energy histogram\n");
  printf("-E <filename>         - Writes energy timeline\n");
  printf("-s <filename>         - Writes spatial end-to-end distribution pdb file\n");
  printf("-r <filename>         - Writes radial end-to-end distribution histogram\n");
  printf("-R <filename>         - Writes radial end-to-end distribution timeline\n");
  printf("-g <filename>         - Writes radius of gyration histogram\n");
  printf("-G <filename>         - Writes radius of gyration timeline\n");
  printf("-S <filename.pdb>     - Writes spatial distribution pdb\n\n");
}

void parsecommandline(int argc, char **argv) {

  strcpy(eh, "");
  strcpy(et, "");
  strcpy(sp, "");
  strcpy(rh, "");
  strcpy(rt, "");
  strcpy(gh, "");
  strcpy(gt, "");

  if (argc < 3) {
    printhelp();
    exit(0);
  }

  strcpy(filename, argv[1]);
  strcpy(source, argv[2]);

  for (int i = 3; i < argc; i++) {

    if (strcmp(argv[i], "-g") == 0) {
      strcpy(gh, argv[i+1]);
      igh = 1;
      i++;
    } else
      if (strcmp(argv[i], "-e") == 0) {
	strcpy(eh, argv[i+1]);
	ieh = 1;
	i++;
      } else 
	if (strcmp(argv[i], "-G") == 0) {
	  strcpy(gt, argv[i+1]);
	  igt = 1;
	  i++;
	} else
	  if (strcmp(argv[i], "-E") == 0) {
	    strcpy(et, argv[i+1]);
	    iet = 1;
	    i++;
	  } else
	    if (strcmp(argv[i], "-r") == 0) {
	      strcpy(rh, argv[i+1]);
	      irh = 1;
	      i++;
	    } else
	      if (strcmp(argv[i], "-R") == 0) {
		strcpy(rt, argv[i+1]);
		irt = 1;
		i++;
	      } else
		if (strcmp(argv[i], "-i") == 0) {
		  sscanf(argv[i+1], "%d", &start);
		  i++;
		} else
		  if (strcmp(argv[i], "-f") == 0) {
		    sscanf(argv[i+1], "%d", &end);
		    i++;
		  } else
      {
	printf("Garbage encountered on command line: %s\n", argv[i]);
	exit(0);
      }

  }
}
