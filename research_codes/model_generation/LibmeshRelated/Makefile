# The location of the mesh library
LIBMESH_DIR = /xlwork6/hgao/IBAMR2017/libmesh-1.0.0/opt

# include the library options determined by configure.  This will
# set the variables INCLUDE and LIBS that we will need to build and
# link with the library.
include $(LIBMESH_DIR)/Make.common


###############################################################################
# File management.  This is where the source, header, and object files are
# defined

#
# source files
srcfiles 	:= $(wildcard *.C)

#
# object files
objects		:= $(patsubst %.C, %.$(obj-suffix), $(srcfiles))
###############################################################################



.PHONY: clean clobber distclean

###############################################################################
# Target:
#
target 	   := ./fibreGen-$(METHOD)


all:: $(target)

# Production rules:  how to make the target - depends on library configuration
$(target): $(objects)
	@echo "Linking "$@"..."
	@$(libmesh_CXX) $(libmesh_CPPFLAGS) $(libmesh_CXXFLAGS) $(objects) -o $@ $(libmesh_LIBS) $(libmesh_LDFLAGS) -L./ -ltet

# Useful rules.
clean:
	@rm -f $(objects) *~

clobber:
	@$(MAKE) clean
	@rm -f $(target) out.gmv

distclean:
	@$(MAKE) clobber
	@rm -f *.o *.g.o *.pg.o .depend

run: $(target)
	@echo "***************************************************************"
	@echo "* Running Example " $(LIBMESH_RUN) $(target) $(LIBMESH_OPTIONS)
	@echo "***************************************************************"
	@echo " "
	@$(LIBMESH_RUN) $(target) $(LIBMESH_OPTIONS)
	@echo " "
	@echo "***************************************************************"
	@echo "* Done Running Example " $(LIBMESH_RUN) $(target) $(LIBMESH_OPTIONS)
	@echo "***************************************************************"


# include the dependency list
include .depend


#
# Dependencies
#
.depend: $(srcfiles) $(LIBMESH_DIR)/include/*/*.h
	@$(perl) $(LIBMESH_DIR)/contrib/bin/make_dependencies.pl -I. $(foreach i, $(wildcard $(LIBMESH_DIR)/include/*), -I$(i)) "-S\$$(obj-suffix)" $(srcfiles) > .depend

###############################################################################
