FC = f77
CC = cc
FFLAGS = -fno-second-underscore
CFLAGS = -Iinclude -Wall

SRC = src
BUILD = build

lib: $(BUILD)/libfilter.a
$(BUILD)/libfilter.a: $(BUILD)/finit.o $(BUILD)/cfkaip.o $(BUILD)/diag.o \
                      $(BUILD)/drnd.o $(BUILD)/lao.o $(BUILD)/pre.o \
                      $(BUILD)/chkinput.o $(BUILD)/psiio.o \
                      $(BUILD)/reai.o $(BUILD)/set_robin_input.o \
                      $(BUILD)/ntpsi.o $(BUILD)/norm_psi.o
	ar rv $(BUILD)/libfilter.a $(BUILD)/*.o

# C files
$(BUILD)/finit.o: $(SRC)/finit.c include/hfile.h include/const.h
	$(CC) $(CFLAGS) -c $(SRC)/finit.c -o $(BUILD)/finit.o

$(BUILD)/pre.o: $(SRC)/pre.c include/const.h
	$(CC) $(CFLAGS) -c $(SRC)/pre.c -o $(BUILD)/pre.o

$(BUILD)/chkinput.o: $(SRC)/chkinput.c include/const.h
	$(CC) $(CFLAGS) -c $(SRC)/chkinput.c -o $(BUILD)/chkinput.o

$(BUILD)/cfkaip.o: $(SRC)/cfkaip.c include/hfile.h
	$(CC) $(CFLAGS) -c $(SRC)/cfkaip.c -o $(BUILD)/cfkaip.o

$(BUILD)/diag.o: $(SRC)/diag.c include/hfile.h
	$(CC) $(CFLAGS) -c $(SRC)/diag.c -o $(BUILD)/diag.o

$(BUILD)/drnd.o: $(SRC)/drnd.c include/hfile.h
	$(CC) $(CFLAGS) -c $(SRC)/drnd.c -o $(BUILD)/drnd.o

$(BUILD)/diag.o: $(SRC)/diag.c include/hfile.h
	$(CC) $(CFLAGS) -c $(SRC)/diag.c -o $(BUILD)/diag.o

$(BUILD)/lao.o: $(SRC)/lao.c include/lao.h
	$(CC) $(CFLAGS) -c $(SRC)/lao.c -o $(BUILD)/lao.o

$(BUILD)/norm_psi.o: $(SRC)/norm_psi.c include/lao.h
	$(CC) $(CFLAGS) -c $(SRC)/norm_psi.c -o $(BUILD)/norm_psi.o

$(BUILD)/ntpsi.o: $(SRC)/ntpsi.c
	$(CC) $(CFLAGS) -c $(SRC)/ntpsi.c -o $(BUILD)/ntpsi.o

$(BUILD)/psiio.o: $(SRC)/psiio.c
	$(CC) $(CFLAGS) -c $(SRC)/psiio.c -o $(BUILD)/psiio.o

# Fortran files
$(BUILD)/reai.o: $(SRC)/reai.F include/const.h
	$(FC) $(FFLAGS) -Iinclude -c $(SRC)/reai.F -o $(BUILD)/reai.o

$(BUILD)/set_robin_input.o: $(SRC)/set_robin_input.F \
    include/filename.h include/spec.h include/logicals.h \
    include/dimensions.h include/f_minmax.h include/gauss.h \
    include/cutprec.h
	$(FC) $(FFLAGS) -Iinclude -c $(SRC)/set_robin_input.F -o $(BUILD)/set_robin_input.o

clean:
	rm -f $(BUILD)/*.o $(BUILD)/libfilter.a

