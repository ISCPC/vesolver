
gcc -std=c99 -O3 -g -fPIC -I../include -I/opt/nec/ve/veos/include -c ../src/vesolver.c
gcc -shared -o libvesolver_api.so  vesolver.o -L/opt/nec/ve/veos/lib64 -Wl,-rpath=/opt/nec/ve/veos/lib64 -lveo
# gcc -o libvesolver_api.so vesolver.o -L/opt/nec/ve/veos/lib64 -Wl,-rpath=/opt/nec/ve/veos/lib64 -lveo
