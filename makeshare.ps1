gcc -Wall -Wextra -pedantic -c -o objects/linalg.o libraries/linalg.c
gcc -Wall -Wextra -pedantic -c -o objects/cycles.o libraries/cycles.c
gcc -Wall -Wextra -pedantic -c -o objects/factors.o libraries/factors.c

gcc -fPIC -shared -o objects/orbitvis.so libraries/orbitvis.c objects/linalg.o objects/cycles.o objects/factors.o