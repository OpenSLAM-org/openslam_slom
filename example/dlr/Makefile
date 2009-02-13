include ../../Makefile.conf

LD = $(CXX)

SRC = ../../src

CPPFLAGS += -I$(SRC)


HEADER = $(SRC)/Estimator.h $(SRC)/types/*.h $(SRC)/tools/*.h $(SRC)/manifolds/*.h

all: dlr_data_set



*.o: $(HEADER) ../../Makefile.conf


dlr_data_set: dlr_data_set.o $(SRC)/Estimator.o
	$(LD) $(LDFLAGS) $^ -o $@ $(LIBS)
