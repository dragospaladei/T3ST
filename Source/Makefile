FC = ifx
FCFLAGS = -qopenmp -O3 -xhost

SOURCES = Import_sims.f90 Constants.f90 rand_omp.f90 Random_number_generation.f90 Drifts_kernels.f90 \
Collisions.f90 Circular_Equilibrium.f90 Dispersion.f90 EFIT_Equilibrium.f90 Equilibrium_magnetic.f90 Magn_Equil_Choice.f90 \
Particle_initialization.f90 Print_messages.f90 Solovev_equilibrium.f90 testing_T3ST.f90 T3ST_main.f90

EXEC = a.out
ONEAPI_ENV_SCRIPT = $(HOME)/intel/oneapi/compiler/2024.1/env/vars.sh

all: $(EXEC)

$(EXEC): $(SOURCES)
	@bash -lc 'source "$(ONEAPI_ENV_SCRIPT)" && \
	           ulimit -s unlimited && \
	           $(FC) $(FCFLAGS) $(SOURCES) -o $(EXEC)'

clean:
	@rm -f $(EXEC)

.PHONY: all clean

