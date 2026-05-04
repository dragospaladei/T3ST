FC = ifx
FCFLAGS = -qopenmp -qopenmp-simd -O3 -xhost
ONEAPI_COMPILER_LIB = $(HOME)/intel/oneapi/compiler/2024.1/lib
LDFLAGS = -Xlinker --disable-new-dtags -Xlinker -rpath -Xlinker $(ONEAPI_COMPILER_LIB)

BUILD_DIR = ../../build
EXEC = $(BUILD_DIR)/T3ST
LEGACY_EXEC = ../../scripts/a.out

ONEAPI_ENV_SCRIPT = $(HOME)/intel/oneapi/compiler/2024.1/env/vars.sh

SOURCES = Import_sims.f90 Constants.f90 Export_files.f90 rand_omp.f90 Random_number_generation.f90 Drifts_kernels.f90 \
Collisions.f90 Circular_Equilibrium.f90 Dispersion.f90 EFIT_Equilibrium.f90 Equilibrium_magnetic.f90 Magn_Equil_Choice.f90 \
Particle_initialization.f90 Print_messages.f90 Solovev_equilibrium.f90 testing_T3ST.f90 T3ST_main.f90

all: $(EXEC)

$(EXEC): $(SOURCES)
	@mkdir -p "$(BUILD_DIR)"
	@bash -lc 'source "$(ONEAPI_ENV_SCRIPT)" && \
	           ulimit -s unlimited && \
	           set -e && \
	           for src in $(SOURCES); do \
	             obj="$(BUILD_DIR)/$${src%.f90}.o"; \
	             $(FC) $(FCFLAGS) -module "$(BUILD_DIR)" -I"$(BUILD_DIR)" -c "$$src" -o "$$obj"; \
	           done && \
	           $(FC) $(FCFLAGS) "$(BUILD_DIR)"/*.o $(LDFLAGS) -o "$(EXEC)"'

clean:
	@rm -f "$(EXEC)" "$(LEGACY_EXEC)" "$(BUILD_DIR)"/*.o "$(BUILD_DIR)"/*.mod

.PHONY: all clean
