CFLAGS   += -Wall -g
CXXFLAGS += $(CFLAGS) 
LDFLAGS  += -lm $(MUSE_LD_FLAGS)

CODELIB = libvader.a

CODEOBJS = test.o advanceBE.o advanceCN.o applyBC.o checkpoint.o driver.o getNextIterate.o init.o ppmExtrap.o rotCurveSpline.o testingMode.o userFunc.o vader.o

AR = ar ruv
RANLIB = ranlib
RM = rm

PROB = none

all: $(CODELIB) 


clean:
	$(RM) -f *.o *.a

userFunc.c: prob/userFunc_$(PROB).c
	cp $< $@

$(CODELIB): $(CODEOBJS)
	$(RM) -f userFunc.c
	$(RM) -f $@
	$(AR) $@ $(CODEOBJS)
	$(RANLIB) $@


.cc.o: $<
	$(CXX) $(CXXFLAGS) -c -o $@ $< 
