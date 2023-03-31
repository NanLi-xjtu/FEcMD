SRCS = mlip.cpp heatdiffusion.cpp singlestep.cpp init.cpp namelist.c utilityfunctions.c measurements.c forces.c heattransport.c vscale.c structure.c eamforce.c eam_alloy.c integral.c nebrList.c tipMaker.c npt.c boundary.c interface.c snap.c

OBJS1 = $(SRCS:.cpp=.o)

OBJS = $(OBJS1:.c=.o)

MLIP_PATH=../lib/mlip-2-master

.PHONY: all clean

all: md

md: $(OBJS) main.o
	g++ -ggdb3 -fopenmp $^ -L$(MLIP_PATH)/lib -l_mlip_interface -L/usr/include/X11/lib -lgfortran -o $@

md_lib: $(OBJS)
	@ar rcs lib/libmd.a $(OBJS)

%.o: %.cpp
	g++ -ggdb3 -fopenmp -c $< -o $@

%.o: %.c
	gcc -ggdb3 -fopenmp -c $< -o $@ -lm

clean:
	rm -f *.o md lib/libmd.a
