#include <stdlib.h>

/* A replacement for system() that uses vfork. S4 seems to get
   confused after using the regular system() but not with this one. */

int vsystem(const char *command)
{
  int pid, w, stat;
  if ((pid = vfork()) == 0) {
    execlp("/bin/posix/sh", "sh", "-c", command, 0);
    _exit(127);
  }
  while ((w = wait(&stat)) != pid && w != -1)
    ;
  if (w == -1)
    stat = -1;
  return stat;
}
