
#ifndef _FILE_CPP
#define _FILE_CPP

#include <stdio.h>
#include <string.h>
#include <unistd.h>

FILE *openfwrite(char *filename) {
  FILE *temp = fopen(filename, "r");
  if (temp != NULL) {
    fclose(temp);
    char s[256];
    strcpy(s, "#");
    strcat(s, filename);
    strcat(s, "#");
    printf("Backing up %s to %s\n", filename, s);
    rename(filename, s);
  }
  return fopen(filename, "w");
}

#endif
