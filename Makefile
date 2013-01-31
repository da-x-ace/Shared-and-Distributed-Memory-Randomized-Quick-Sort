all:
	cilk++ -m64 -fPIC -shared -o libac.so 1c.cilk
	mpicxx 2b.cpp -Wl,-rpath=. -L. -lac -o prob
