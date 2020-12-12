#include <stdio.h>
#include <stdlib.h>

#include "bdna/bdna.h"

main() {
  bdna b(146);
  b.readconfig("out.steps");
  printf("Twist = %lf, Writhe = %lf\n", b.calculate_twist_closed(), b.calculate_writhe_closed());
}
