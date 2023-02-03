FC = gfortran
FLAGS = -O3  # -fbounds-check
TARGET = test 

SOURCES = $(wildcard *.f90)
OBJS = $(SOURCES:.f90=.o)

%.o: %.f90 
	$(FC) $(FLAGS) -c  $<

all: $(OBJS)
	$(FC) -o $(TARGET) $(OBJS) $(LIBS)

clean:
	rm *.o *.mod $(TARGET)

include make.dep
