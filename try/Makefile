CC = f77
CFLAGS = -Wall -O2 -g
LIB = ../lib/build/libfilter.a
SRC = src

.PHONY: all clean

all: test_fd_base test_lanczos_first test_flset

TEST_UTIL_OBJS = $(SRC)/fd_test_utils.o

TEST_LANCZ_OBJS = $(SRC)/fd_test_lanczos.o

# Object files
TEST_UTIL_OBJS = $(SRC)/fd_test_utils.o

$(SRC)/fd_test_utils.o: $(SRC)/fd_test_utils.c $(SRC)/fd_test_utils.h
	$(CC) $(CFLAGS) -c $< -o $@

$(SRC)/fd_test_lanczos.o: $(SRC)/fd_test_lanczos.c $(SRC)/fd_test_lanczos.h $(SRC)/fd_test_utils.h
	$(CC) $(CFLAGS) -c $< -o $@

test_fd_base.o: $(SRC)/test_fd_base.c $(SRC)/fd_test_utils.h $(SRC)/fd_test_lanczos.h
	$(CC) $(CFLAGS) -c $< -o $@

test_lanczos_first.o: $(SRC)/test_lanczos_first.c $(SRC)/fd_test_utils.h $(SRC)/fd_test_lanczos.h
	$(CC) $(CFLAGS) -c $< -o $@

test_flset.o: $(SRC)/test_flset.c $(SRC)/fd_test_utils.h $(SRC)/fd_test_lanczos.h
	$(CC) $(CFLAGS) -c $< -o $@

test_fd_base: test_fd_base.o $(TEST_UTIL_OBJS) $(LIB) $(TEST_LANCZ_OBJS)
	$(CC) $(CFLAGS) -o $@ test_fd_base.o $(TEST_UTIL_OBJS) $(TEST_LANCZ_OBJS) $(LIB) -llapack -lblas -lgfortran -lm

test_lanczos_first: test_lanczos_first.o $(TEST_UTIL_OBJS) $(TEST_LANCZ_OBJS) $(LIB)
	$(CC) $(CFLAGS) -o $@ test_lanczos_first.o $(TEST_UTIL_OBJS) $(TEST_LANCZ_OBJS) $(LIB) -llapack -lblas -lgfortran -lm

test_flset: test_flset.o $(TEST_UTIL_OBJS) $(TEST_LANCZ_OBJS) $(LIB)
	$(CC) $(CFLAGS) -o $@ test_flset.o $(TEST_UTIL_OBJS) $(TEST_LANCZ_OBJS) $(LIB) -llapack -lblas -lgfortran -lm

# Build library if needed
$(LIB):
	$(MAKE) -C ../lib

# Clean
clean:
	rm -f $(SRC)/*.o test_fd_base test_lanczos_first test_flset

