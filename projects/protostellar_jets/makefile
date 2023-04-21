ifndef AMRVAC_DIR
$(error AMRVAC_DIR is not set)
endif

ARCH = openmpi_run
ifdef ARCH
	ARCH := $(ARCH)
else
  ifndef ARCH_LOC
	ARCH := default
  else
	ARCH := $(ARCH_LOC)
  endif
endif

ifdef ndim
	NDIM := $(ndim)
else
  ifndef NDIM_LOC
	NDIM := 2
  else
	NDIM := $(NDIM_LOC)
  endif
endif

# By exporting these can be used when building libamrvac
export ARCH NDIM

OBJ_DIR := Objects
SRC_DIR := $(AMRVAC_DIR)/src
LIB_DIR := $(AMRVAC_DIR)/lib/$(NDIM)d_$(ARCH)
LIB_MAKE := $(AMRVAC_DIR)/arch/lib.make
LIB_AMRVAC := $(LIB_DIR)/libamrvac.a
LIB_USER := lib_user
SRC_DIR += $(OBJ_DIR)

# These are used for compilation
INC_DIRS := $(LIB_DIR)
LIB_DIRS := $(LIB_DIR)
LIBS := amrvac



.PHONY: all clean allclean force hdr



all: | hdr amrvac

# copy amrvac.h (in order to use the std preprocessor in the main code files, e.g. twofl); create the file if it does not exist
hdr:
ifeq ("$(wildcard amrvac.h)","")
	touch amrvac.h
endif
	@mkdir -p $(LIB_DIR)	# Prevent error message
	rsync -c amrvac.h $(LIB_DIR)/amrvac.h


# Include architecture and compilation rules
include $(AMRVAC_DIR)/arch/$(ARCH).defs
include $(AMRVAC_DIR)/arch/rules.make

# Optionally include a local user makefile
-include local.make

# Where to find amrvac.t
vpath %.t $(SRC_DIR)

# Keep mod_usr.f for inspection
.PRECIOUS: mod_usr.f $(OBJ_DIR)/mod_obj_global_parameters.f \
	$(OBJ_DIR)/mod_obj_types.f \
	$(OBJ_DIR)/mod_obj_read_parameters.f \
	$(OBJ_DIR)/mod_obj_input_output.f

# Intermediate files are removed
.INTERMEDIATE: amrvac.o mod_usr.o mod_usr.mod

# Always try to build/update the amrvac library
$(LIB_AMRVAC): force
	@mkdir -p $(LIB_DIR)
	@mkdir -p $(LIB_USER)
	$(MAKE) -C $(LIB_DIR) -f $(LIB_MAKE)

clean:
	@echo $(SRC_DIR)
	@echo 'Cleaning local objects ("make allclean" cleans libamrvac)'
	@mkdir -p $(LIB_DIR)	# Prevent error message
	$(RM) amrvac amrvac.o mod_usr.o mod_usr.mod\
		mod_obj_**.o mod_obj_*.mod\
	    $(OBJ_DIR)/*.o     \
		$(OBJ_DIR)/*.mod   \
		$(OBJ_DIR)/*.f *.f


# Also clean the library
allclean: clean
	@echo 'Cleaning libamrvac'
	@mkdir -p $(LIB_DIR)	# Prevent error message
	$(MAKE) -C $(LIB_DIR) -f $(LIB_MAKE) clean

# Dependencies
amrvac: mod_usr.o amrvac.o $(OBJ_DIR)/mod_obj_global_parameters.o \
	$(OBJ_DIR)/mod_obj_types.o \
	$(OBJ_DIR)/mod_obj_read_parameters.o \
	$(OBJ_DIR)/mod_obj_input_output.o


amrvac.o mod_usr.o: $(LIB_AMRVAC)
amrvac.o: mod_usr.o

mod_obj_global_parameters.o : $(LIB_AMRVAC)

mod_obj_types.o : $(LIB_AMRVAC) \
	$(OBJ_DIR)/mod_obj_global_parameters.o

mod_obj_read_parameters.o : $(LIB_AMRVAC) \
	$(OBJ_DIR)/mod_obj_global_parameters.o \
	$(OBJ_DIR)/mod_obj_types.o

mod_obj_input_output.o : $(LIB_AMRVAC) \
	$(OBJ_DIR)/mod_obj_global_parameters.o \
	$(OBJ_DIR)/mod_obj_types.o \
	$(OBJ_DIR)/mod_obj_read_parameters.o

mod_usr.o: mod_usr.f amrvac.h \
	$(OBJ_DIR)/mod_obj_global_parameters.o \
	$(OBJ_DIR)/mod_obj_types.o \
	$(OBJ_DIR)/mod_obj_read_parameters.o \
	$(OBJ_DIR)/mod_obj_input_output.o
