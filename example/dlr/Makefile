include ../Make.include


HEADER += dlr_parser.h dlr_types.h ../tools.h

all: dlr_data_set dlr_data_set_calib



*.o: $(HEADER) ../../Makefile.conf


dlr_data_set_calib.o: dlr_data_set.cpp $(HEADER) ../../Makefile.conf
	$(CXX) $(CPPFLAGS) -DDLR_CALIBRATE -c -o $@  $<

dlr_data_set: dlr_data_set.o $(SRC)/Estimator.o
	$(LD) $(LDFLAGS) $^ -o $@ $(LIBS)

dlr_data_set_calib: dlr_data_set_calib.o $(SRC)/Estimator.o
	$(LD) $(LDFLAGS) $^ -o $@ $(LIBS)
	