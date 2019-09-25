# standard amuse configuration include
# config.mk will be made after ./configure has run
AMUSE_DIR?=../../amuse
-include $(AMUSE_DIR)/config.mk

MPICXX   ?= mpicxx

CFLAGS   += -Wall -g
CXXFLAGS += $(CFLAGS) 
LDFLAGS  += -lm -lgsl -lgslcblas $(MUSE_LD_FLAGS)

OBJS = interface.o

CODELIB = src/libvader.a

CODE_GENERATOR = $(AMUSE_DIR)/build.py

PROB = none

all: vader_worker 

clean:
	$(RM) -f *.so *.o *.pyc worker_code.cc worker_code.h 
	$(RM) *~ vader_worker worker_code.cc
	make -C src clean

$(CODELIB):
	make -C src all PROB=$(PROB)

worker_code.cc: interface.py
	$(CODE_GENERATOR) --type=c interface.py vaderInterface -o $@

worker_code.h: interface.py
	$(CODE_GENERATOR) --type=H interface.py vaderInterface -o $@

vader_worker: worker_code.cc worker_code.h $(CODELIB) $(OBJS)
	$(MPICXX) $(CXXFLAGS) $< $(OBJS) $(CODELIB) -o $@ $(LDFLAGS)

.cc.o: $<
	$(CXX) $(CXXFLAGS) -c -o $@ $< 
