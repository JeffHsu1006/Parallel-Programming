NVFLAGS  := -std=c++11 -O3 -Xptxas="-v" -arch=sm_61 
LDFLAGS  := -lm
NVCCLDFLAGS  = -L/home/pp19/share/lib -Xlinker=-rpath,/home/pp19/share/lib -llab5

EXES     := seq hw4-1

alls: $(EXES)

clean:
	rm -f $(EXES)

seq: seq.cc
	g++ $(CXXFLAGS) -o $@ $?

hw4-1: hw4-1.cu
	nvcc $(NVFLAGS) $(NVCCLDFLAGS) -o $@ $?
