MDIR := $(shell pwd)

CC = gcc
LIB_LINK = -lgsl -lgslcblas -lm -lfftw3
FLAGS = -fPIC -Wall -O3 -ffast-math
BUILD_DIR = $(MDIR)/build
LIB = libfftlogx
SRC_DIR = $(MDIR)/src
_SRC = \
	cfftlog.c \
	utils.c \
	utils_complex.c \

SRC = $(addprefix $(SRC_DIR)/,$(_SRC))


all: $(LIB).so

$(LIB).so:
	if ! [ -e $(BUILD_DIR) ]; then mkdir $(BUILD_DIR); fi;
	$(CC) -shared -o $(BUILD_DIR)/$(LIB).so $(LIB_LINK) $(FLAGS) $(SRC)

clean:
	rm $(BUILD_DIR)/*.so