cdir=${HOME}/clib/
CFLAGS=-Wall -W \
        -Wmissing-prototypes -Wstrict-prototypes \
        -Wconversion -Wshadow \
        -Wpointer-arith -Wcast-qual -Wcast-align \
        -Wwrite-strings -Wnested-externs \
        -fno-common -Dinline= -g -O4 -I/usr/local/include/ \

#CFLAGS=-O4 -DHAVE_INLINE -DGSL_RANGE_CHECK_OFF -I/usr/local/include/

cdirfiles=${cdir}matrix.o ${cdir}ran.o ${cdir}distr.o ${cdir}ps.o ${cdir}array.o ${cdir}bfgs.o

files=main.o VD.o 611wrapper.o 611new.o post-opt.o mvtdstpack.o UO.o ${cdirfiles}


all: ${files}
	gcc ${CFLAGS} ${files} -I${cdir}  -L/usr/local/lib/ -L/opt/local/lib/ -lgsl -lgslcblas -lf2c -lm -o run


main.o : main.c
	gcc ${CFLAGS} -c main.c -I${cdir} 
VD.o : VD.c VD.h
	gcc ${CFLAGS} -c VD.c -I${cdir} 
post-opt.o : post-opt.c post-opt.h
	gcc ${CFLAGS} -c post-opt.c -I${cdir} 


611new.o : 611new.c 611new.h
	gcc ${CFLAGS} -c 611new.c -I${cdir} 
611wrapper.o : 611wrapper.c 611wrapper.h
	gcc ${CFLAGS} -c 611wrapper.c -I${cdir} 
mvtdstpack.o : mvtdstpack.c mvtdstpack.h
	gcc ${CFLAGS} -c mvtdstpack.c -I${cdir} 
UO.o : UO.c UO.h
	gcc ${CFLAGS} -c UO.c -I${cdir} 





${cdirfiles}: 
	make -C ${cdir} CFLAGS='${CFLAGS}'

clean:
	rm -rf *.o run ${cdir}*.o sim* *.out cum*
