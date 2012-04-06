#include <dl.h>

#define RTLD_LAZY (BIND_DEFERRED | BIND_NONFATAL)

void *dlopen(const char *, int);
void *dlsym(void *, const char *);
int dlclose(void *);
char *dlerror(void);

