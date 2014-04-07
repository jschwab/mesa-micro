#################################################################

# STEP 1: get the mesa makefile_header

include $(MESA_DIR)/utils/makefile_header

#################################################################

# STEP 2: build the test

MICRO = micro
MICRO_OBJS =  micro_support.o micro.o

all: $(MICRO)

$(MICRO) : $(MICRO_OBJS)
	$(LOADER) $(FCopenmp) -o $(MICRO) $(MICRO_OBJS) \
	-L$(MESA_LIB_DIR) $(LOAD_MESA_MICRO)

# for more options, see utils/makefile_header

#################################################################

MY_FC_FLAGS = $(FCfree)

%.o: %.f
	$(FC) $(FCchecks) $(MY_FC_FLAGS) -I$(MESA_INCLUDE_DIR) -c $<

clean:
	-@rm -f *.o *.mod $(MICRO)
