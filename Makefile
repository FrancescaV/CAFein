OBJS = dma.o SteffenInterp.o gsl_odeiv_evolve_Riccati.o matrix_operations.o RiccatiMatrices_operations.o \
eigenfunction_operations.o IOfiles.o initialConditions.o readInputParameters.o TidalParameters_operations.o \
numericalIntegration.o
HEADERS = gsl_odeiv_evolve_Riccati.h dma.h SteffenInterp.h params_integrator_V_nonAd.h params_integrator_adiabatic.h \
 matrix_operations.h RiccatiMatrices_operations.h eigenfunction_operations.h IOfiles.h \
 initialConditions.h readInputParameters.h TidalParameters_operations.h numericalIntegration.h \
 params_integrator_R_nonAd.h rkf45_state.h params_integrator_V_nonAd_rkf45.h
CC = g++
DEBUG = -g
CFLAGS = -Wall -Wunused -O3 -c $(DEBUG)
LFLAGS = -Wall -Wunused -O3 $(DEBUG)
EXECUTABLE=CAFein
LIBS=-lgsl -lgslcblas -lm

$(EXECUTABLE): $(EXECUTABLE).o $(OBJS)
	$(CC) $(LFLAGS) $(EXECUTABLE).o $(OBJS) $(LIBS) -o $(EXECUTABLE)
	
$(EXECUTABLE).o: $(EXECUTABLE).c $(HEADERS)
	$(CC) $(CFLAGS) $(EXECUTABLE).c

dma.o: dma.h dma.c
	$(CC) $(CFLAGS) dma.c

SteffenInterp.o: SteffenInterp.h SteffenInterp.c dma.h
	$(CC) $(CFLAGS) SteffenInterp.c

gsl_odeiv_evolve_Riccati.o: gsl_odeiv_evolve_Riccati.h RiccatiMatrices_operations.h \
params_integrator_V_nonAd.h params_integrator_adiabatic.h params_integrator_R_nonAd.h SteffenInterp.h \
eigenfunction_operations.h matrix_operations.h IOfiles.h params_integrator_V_nonAd_rkf45.h\
gsl_odeiv_evolve_Riccati.c 
	$(CC) $(CFLAGS) gsl_odeiv_evolve_Riccati.c

matrix_operations.o: matrix_operations.h IOfiles.h matrix_operations.c 
	$(CC) $(CFLAGS) matrix_operations.c

RiccatiMatrices_operations.o: RiccatiMatrices_operations.h matrix_operations.h IOfiles.h RiccatiMatrices_operations.c 
	$(CC) $(CFLAGS) RiccatiMatrices_operations.c

eigenfunction_operations.o: eigenfunction_operations.h SteffenInterp.h IOfiles.h dma.h numericalIntegration.h rkf45_state.h\
eigenfunction_operations.c
	$(CC) $(CFLAGS) eigenfunction_operations.c
	
IOfiles.o: IOfiles.h IOfiles.c
	$(CC) $(CFLAGS) IOfiles.c

initialConditions.o: initialConditions.h RiccatiMatrices_operations.h initialConditions.c
	$(CC) $(CFLAGS) initialConditions.c

readInputParameters.o: readInputParameters.h IOfiles.h readInputParameters.c
	$(CC) $(CFLAGS) readInputParameters.c

TidalParameters_operations.o: TidalParameters_operations.h SteffenInterp.h IOfiles.h numericalIntegration.h dma.h TidalParameters_operations.c
	$(CC) $(CFLAGS) TidalParameters_operations.c

numericalIntegration.o: numericalIntegration.h IOfiles.h numericalIntegration.c
	$(CC) $(CFLAGS) numericalIntegration.c

.PHONY = clean

clean:
	\rm *.o $(EXECUTABLE)
