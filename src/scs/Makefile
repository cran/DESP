# MAKEFILE for scs
include scs.mk

OBJECTS_SCS = src/scs.o src/util.o src/cones.o src/cs.o src/linAlg.o src/ctrlc.o src/scs_version.o

SRC_FILES = $(wildcard src/*.c)
INC_FILES = $(wildcard include/*.h)

AMD_SOURCE = $(wildcard $(DIRSRCEXT)/amd_*.c)
DIRECT_OBJECTS_SCS = $(DIRSRCEXT)/ldl.o $(AMD_SOURCE:.c=.o) 
TARGETS = $(OUT)/demo_direct $(OUT)/demo_indirect $(OUT)/demo_SOCP_indirect $(OUT)/demo_SOCP_direct

.PHONY: default 

default: $(TARGETS) $(LIB)/libscsdir.a $(LIB)/libscsindir.a $(LIB)/libscsdir.$(SHARED) $(LIB)/libscsindir.$(SHARED)
	@echo "**********************************************************************************"
	@echo "Successfully compiled scs, copyright Brendan O'Donoghue 2012."
	@echo "To test, type '$(OUT)/demo_direct' or '$(OUT)/demo_indirect'."
	@echo "**********************************************************************************"
ifneq ($(USE_LAPACK), 0)
	@echo "Compiled with blas and lapack, can solve LPs, SOCPs, SDPs, and ECPs"
else
	@echo "NOT compiled with blas/lapack, cannot solve SDPs (can solve LPs, SOCPs, and ECPs)."
	@echo "To solve SDPs, install blas and lapack, then edit scs.mk to set USE_LAPACK=1"
	@echo "and point to the library install locations, and recompile with 'make purge', 'make'."
endif
	@echo "**********************************************************************************"

src/scs.o	: $(SRC_FILES) $(INC_FILES)
src/util.o	: src/util.c include/util.h include/constants.h
src/cones.o	: src/cones.c include/cones.h
src/cs.o	: src/cs.c include/cs.h
src/linAlg.o: src/linAlg.c include/linAlg.h
src/ctrl.o  : src/ctrl.c include/ctrl.h
src/scs_version.o: src/scs_version.c include/constants.h

$(DIRSRC)/private.o: $(DIRSRC)/private.c  $(DIRSRC)/private.h
$(INDIRSRC)/indirect/private.o: $(INDIRSRC)/private.c $(INDIRSRC)/private.h
$(LINSYS)/common.o: $(LINSYS)/common.c $(LINSYS)/common.h

$(LIB)/libscsdir.a: $(OBJECTS_SCS) $(DIRSRC)/private.o $(DIRECT_OBJECTS_SCS) $(LINSYS)/common.o
	mkdir -p $(LIB)
	$(ARCHIVE) $(LIB)/libscsdir.a $^
	- $(RANLIB) $(LIB)/libscsdir.a

$(LIB)/libscsindir.a: $(OBJECTS_SCS) $(INDIRSRC)/private.o $(LINSYS)/common.o
	mkdir -p $(LIB)
	$(ARCHIVE) $(LIB)/libscsindir.a $^
	- $(RANLIB) $(LIB)/libscsindir.a

$(LIB)/libscsdir.$(SHARED): $(OBJECTS_SCS) $(DIRSRC)/private.o $(DIRECT_OBJECTS_SCS) $(LINSYS)/common.o
	mkdir -p $(LIB)
	$(CC) -shared -o $@ $^ $(LDFLAGS)

$(LIB)/libscsindir.$(SHARED): $(OBJECTS_SCS) $(INDIRSRC)/private.o $(LINSYS)/common.o
	mkdir -p $(LIB)
	$(CC) -shared -o $@ $^ $(LDFLAGS)

$(OUT)/demo_direct: examples/c/demo.c $(LIB)/libscsdir.a
	mkdir -p $(LIB)
	mkdir -p $(OUT)
	$(CC) $(CFLAGS) -DDEMO_PATH="\"$(CURDIR)/examples/raw/demo_data\"" $^ -o $@ $(LDFLAGS)

$(OUT)/demo_indirect: examples/c/demo.c $(LIB)/libscsindir.a
	mkdir -p $(LIB)
	mkdir -p $(OUT)
	$(CC) $(CFLAGS) -DDEMO_PATH="\"$(CURDIR)/examples/raw/demo_data\"" $^  -o $@ $(LDFLAGS)

$(OUT)/demo_SOCP_direct: examples/c/randomSOCPProb.c $(LIB)/libscsdir.$(SHARED)
	mkdir -p $(LIB)
	mkdir -p $(OUT)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

$(OUT)/demo_SOCP_indirect: examples/c/randomSOCPProb.c $(LIB)/libscsindir.$(SHARED)
	mkdir -p $(LIB)
	mkdir -p $(OUT)
	$(CC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

.PHONY: clean purge
clean:
	@rm -rf $(TARGETS) $(OBJECTS_SCS) $(DIRECT_OBJECTS_SCS) $(LINSYS)/common.o $(DIRSRC)/private.o $(INDIRSRC)/private.o
	@rm -rf $(OUT)/*.dSYM
	@rm -rf .idea
purge: clean 
	@rm -rf $(LIB)
	@rm -rf $(OUT)

