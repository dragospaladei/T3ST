FC = ifx
FCFLAGS = -qopenmp -qopenmp-simd -O3 -xhost
FCFLAGS_PARTICLE_INITIAL_CONDITIONS = -qopenmp -qopenmp-simd -O1 -xhost
ONEAPI_COMPILER_LIB = $(HOME)/intel/oneapi/compiler/2024.1/lib
LDFLAGS = -Xlinker --disable-new-dtags -Xlinker -rpath -Xlinker $(ONEAPI_COMPILER_LIB)

BUILD_DIR = ../../build
EXEC = $(BUILD_DIR)/T3ST
LEGACY_EXEC = ../../scripts/a.out

ONEAPI_ENV_SCRIPT = $(HOME)/intel/oneapi/compiler/2024.1/env/vars.sh

SOURCES = Import_sims.f90 Run_parameters.f90 Export_files.f90 Thread_rng.f90 Turbulence_random_fields.f90 Gyro_center_kernels.f90 \
Circular_Equilibrium.f90 Dispersion.f90 EFIT_Equilibrium.f90 Equilibrium_magnetic.f90 Equilibrium_dispatch.f90 \
Particle_initial_conditions.f90 Warnings_printing.f90 Solovev_equilibrium.f90 testing_T3ST.f90 T3ST_main.f90

OBJECTS = $(addprefix $(BUILD_DIR)/,$(SOURCES:.f90=.o))

all: $(EXEC)

$(EXEC): $(SOURCES)
	@mkdir -p "$(BUILD_DIR)"
	@bash -lc 'source "$(ONEAPI_ENV_SCRIPT)" && \
	           ulimit -s unlimited && \
	           set -e && \
	           for src in $(SOURCES); do \
	             obj="$(BUILD_DIR)/$${src%.f90}.o"; \
	             flags="$(FCFLAGS)"; \
	             if [ "$$src" = "Particle_initial_conditions.f90" ]; then \
	               flags="$(FCFLAGS_PARTICLE_INITIAL_CONDITIONS)"; \
	             fi; \
	             $(FC) $$flags -module "$(BUILD_DIR)" -I"$(BUILD_DIR)" -c "$$src" -o "$$obj"; \
	           done && \
	           $(FC) $(FCFLAGS) $(OBJECTS) $(LDFLAGS) -o "$(EXEC)"'

clean:
	@rm -f "$(EXEC)" "$(LEGACY_EXEC)" "$(BUILD_DIR)"/*.o "$(BUILD_DIR)"/*.mod

.PHONY: all clean
