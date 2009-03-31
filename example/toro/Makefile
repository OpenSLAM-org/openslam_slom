include ../Make.include

HEADER += 

all: relation2d relation3d



*.o: $(HEADER) ../../Makefile.conf

relation2d: relation2d.o $(SRC)/Estimator.o
	$(LD) $(LDFLAGS) $^ -o $@ $(LIBS)

relation3d: relation3d.o $(SRC)/Estimator.o
	$(LD) $(LDFLAGS) $^ -o $@ $(LIBS)
