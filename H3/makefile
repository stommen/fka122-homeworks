#  CC = gcc

#  LOCAL=${HOME}/.local/
#  CFLAGS = \
#  	 -g \
#  	 -pedantic \
#  	 -Wall \
#  	 -Werror \
#  	 -Iinclude \
#  	 -Wl,-rpath,${LOCAL}/lib \
#  	 -L${LOCAL}/lib \
#  	 -I${LOCAL}/include

	 
#  CFLAGS_OPT = \
#  	     -O2 \
#  	     -march=native

#  LIB = \
#        -lm \
#        -lgsl \
#        -lgslcblas

#  OBJ = \
#        obj/tools.o \
#        obj/run.o
#  MAIN = \
#         obj/main.o

#  ifeq ($(MAKECMDGOALS),test)
#  OBJ_TEST = $(patsubst obj/%.o,obj_test/%.o,$(OBJ))
#  -include unit-test/test.mk
#  else
#  CFLAGS += $(CFLAGS_OPT)
#  endif

#  H3: obj _H3

#  _H3: $(MAIN) $(OBJ)
#  	$(CC) $(CFLAGS) $^ -o H3 $(LIB)

#  obj/%.o: src/%.c
#  	$(CC) -MMD -c $(CFLAGS) $< -o $@ 

#  obj:
#  	mkdir -p obj

#  clean:
#  	find -iname "*.o" -exec rm {} \;
#  	find -iname "*.d" -exec rm {} \;
#  	rm -f H3 run-test
#  	rm -rf obj obj_test

#  .PHONY: clean



 #------------------------ JONATANS MAKE ----------------------------------------------
 CC = clang
 LIBDIR = ${LOCAL}/lib
 LOCAL = ${HOME}/.local

 CFLAGS = \
 	-pedantic \
 	-Wall \
 	-Werror \
 	-g \
 	-Iinclude \
 	-MMD \
 	-I${LOCAL}/include \
 	-I/opt/homebrew/include \
 	-Isrc \

 LDFLAGS = \
 	-Wl,-rpath,/opt/homebrew/lib \
 	-L/opt/homebrew/lib

 LIB = \
 	/opt/homebrew/lib/libgsl.dylib \
 	/opt/homebrew/lib/libgslcblas.dylib \
 	-lm

 OBJ = \
       obj/tools.o \
       obj/run.o

#  CXXFLAGS = \
#  	   -DUSE_TABLE_FOR_RADIAL_FUNCTIONS \
#  	   -O3 \
#  	   -fPIC \
#  	   -std=c++11

 MAIN = \
        obj/main.o

 ifeq ($(MAKECMDGOALS),test)
 OBJ_TEST = $(patsubst obj/%.o,obj_test/%.o,$(OBJ))
 -include unit-test/test.mk
 else
 CFLAGS += $(CFLAGS_OPT)
 endif

 H3: obj _H3

 
 _H3: $(MAIN) $(OBJ)	
	$(CC) $(CFLAGS) $^ -o H3 $(LIB)

 obj/%.o: src/%.c
	$(CC) -MMD -c $(CFLAGS) $< -o $@ 

 obj:
	mkdir -p obj

clean:
	find . -name "*.o" -exec rm {} \;
	find . -name "*.d" -exec rm {} \;
	rm -f H3 run-test
	rm -rf obj obj_test ${LIBDIR}

.PHONY: clean