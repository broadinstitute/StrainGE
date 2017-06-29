#include <stdlib.h>
#include <stdio.h>
#include <string.h>

int
main(int argc, char **argv)
{
  int count = 0;
  FILE *out;
  char line[10000];

  if (argc != 4) {
    fprintf(stderr, "usage: splitfastq input.fq.gz out1.fq.gz out2.fq.gz\n");
    exit(1);
  }

  strcpy(line, "zcat <");
  strcat(line, argv[1]);
  FILE *in = popen(line, "r");

  strcpy(line, "gzip >");
  strcat(line, argv[2]);
  FILE *out1 = popen(line, "w");

  strcpy(line, "gzip >");
  strcat(line, argv[3]);
  FILE *out2 = popen(line, "w");
    
  while (fgets(line, sizeof(line), in)) {
    out = (count & 1) ? out2 : out1;
    if (line[0] != '@') {
      fprintf(stderr, "no @: %s", line);
      exit(1);
    }
    fputs(line, out);
    fgets(line, sizeof(line), in);
    fputs(line, out);
    fgets(line, sizeof(line), in);
    if (line[0] != '+') {
      fprintf(stderr, "no +: %s", line);
      exit(1);
    }
    fputs(line, out);
    fgets(line, sizeof(line), in);
    fputs(line, out);
    ++count;
  }
}
