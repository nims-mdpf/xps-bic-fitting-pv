CC = g++
LINK = $(CC)
RM = rm

TARGET = activeshirley_gcc.exe
OBJ = util.obj settings.obj detectpeaks.obj calcinitialbackground.obj minimumpoints.obj optimizeparameters.obj output.obj optimizeflow.obj mainflow.obj main.obj

SRC = $(OBJ:%obj=%.cpp)

# CPPFLAGS = -g
CPPFLAGS = -O3
# CPPFLAGS = -O3 -march=native
CPPFLAGS += -std=c++11 -fopenmp
LFLAGS = -lm
LIBS = -L/usr/lib -llapack -lblas -fopenmp

PACKDIR = ../package

all: $(TARGET)

$(TARGET): $(OBJ)
	$(LINK) -o $@ $(LFLAGS) $(OBJ) $(LIBS)

%.obj: %.cpp
	$(CC) -o $@ -c $(CPPFLAGS) $<

clean:
	$(RM) -f $(TARGET) $(OBJ)

test: $(TARGET)
	cd C1s_CuO; export OMP_NUM_THREADS=16; /usr/bin/time ./test_single.sh
#	cd MIDATA001.107; export OMP_NUM_THREADS=16; /usr/bin/time ./test_fnd.sh
#	cd sample/As2O3_O1s; ../../$(TARGET) As2O3_O1s.csv

package: $(TARGET)
#	rm -rf $(PACKDIR)
	mkdir -p $(PACKDIR)/src/
	cp -a $(TARGET) automatic_xps_peak_separation_single.py make_result_figures_ppt2.py myfunc.py plot_result_shirley_and_linear.py recipe.txt settings.inp exec_single.sh $(PACKDIR)/src

.PHONY: clean
