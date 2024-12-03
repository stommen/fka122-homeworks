CC = gcc
LOCAL=${HOME}/.local

CFLAGS = \
	 -g \
	 -pedantic \
	 -Werror \
	 -Wall \
	 -Iinclude \
	 -I${LOCAL}/include

LDFLAGS += \
    -L${LOCAL}/lib \
    -Llib \
    -Wl,-rpath,${LOCAL}/lib \
    -Wl,-rpath,lib

CFLAGS_OPT = \
	     -O2 \
	     -march=native

LIB = \
      -lm \
      -lgsl \
      -lgslcblas

OBJ = \
      obj/tools.o \
	  obj/lattice.o \
	  obj/run.o

CXXFLAGS = \
	   -DUSE_TABLE_FOR_RADIAL_FUNCTIONS \
	   -O3 \
	   -fPIC \
	   -std=c++11

MAIN = \
       obj/main.o

ifeq ($(MAKECMDGOALS),test)
OBJ_TEST = $(patsubst obj/%.o,obj_test/%.o,$(OBJ))
-include unit-test/test.mk
else
CFLAGS += $(CFLAGS_OPT)
endif

_H2: obj H2

obj/%.o: src/%.c
	$(CC) -MMD -c $(CFLAGS) $< -o $@

H2: $(MAIN) $(OBJ)
	$(CC) $(CFLAGS) $^ $(LIB) -o $@

obj:
	mkdir -p obj
	mkdir -p obj_nep

clean:
	find . -name "*.o" -exec rm {} \;
	find . -name "*.d" -exec rm {} \;
	rm -f H2 run-test
	rm -rf obj obj_test 

.PHONY: clean


# ------------------------ JONATANS MAKE ----------------------------------------------
# CC = clang
# LIBDIR = ${LOCAL}/lib
# LOCAL = ${HOME}/.local

# CFLAGS = \
# 	-pedantic \
# 	-Wall \
# 	-Werror \
# 	-g \
# 	-Iinclude \
# 	-MMD \
# 	-I${LOCAL}/include \
# 	-I/opt/homebrew/include \
# 	-Isrc \

# LDFLAGS = \
# 	-Wl,-rpath,/opt/homebrew/lib \
# 	-L/opt/homebrew/lib

# LIB = \
# 	/opt/homebrew/lib/libgsl.dylib \
# 	/opt/homebrew/lib/libgslcblas.dylib \
# 	-lm

# OBJ = \
#       obj/tools.o \
#       obj/run.o

# # OBJ_LIB = \
# # 	obj_nep/nep.o \
# # 	obj_nep/interface.o

# CXXFLAGS = \
# 	   -DUSE_TABLE_FOR_RADIAL_FUNCTIONS \
# 	   -O3 \
# 	   -fPIC \
# 	   -std=c++11

# MAIN = \
#        obj/main.o

# ifeq ($(MAKECMDGOALS),test)
# OBJ_TEST = $(patsubst obj/%.o,obj_test/%.o,$(OBJ))
# -include unit-test/test.mk
# else
# CFLAGS += $(CFLAGS_OPT)
# endif

# _H2: obj H2

# obj/%.o: src/%.c
# 	$(CC) -MMD -c $(CFLAGS) $< -o $@

# H2: $(MAIN) $(OBJ)
# 	$(CC) $(CFLAGS) $^ $(LIB) -o $@


# obj:
# 	mkdir -p obj
# 	mkdir -p ${LIBDIR}

# clean:
# 	find . -name "*.o" -exec rm {} \;
# 	find . -name "*.d" -exec rm {} \;
# 	rm -f H2 run-test
# 	rm -rf obj obj_test obj_nep ${LIBDIR}

# .PHONY: clean
