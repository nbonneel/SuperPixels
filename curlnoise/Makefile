# This is for GNU make; other versions of make may not run correctly.

MAIN_PROGRAM = curlnoise

# use this for the 2D examples
#SRC = main2.cpp gluvi.cpp noise.cpp example_2dwake.cpp example_2dmouse.cpp example_2dfancy.cpp

# use this for the 3D example
SRC = main3.cpp gluvi.cpp noise.cpp example_3dplume.cpp

include Makefile.defs

# object files
RELEASE_OBJ = $(patsubst %.cpp,obj/%.o,$(notdir $(SRC)))
DEBUG_OBJ = $(patsubst %.cpp,obj_debug/%.o,$(notdir $(SRC)))

# how to make the main target (debug mode, the default)
$(MAIN_PROGRAM): $(DEBUG_OBJ)
	$(LINK) $(DEBUG_LINKFLAGS) -o $@ $^ $(LINK_LIBS)

# how to make the main target (release mode)
$(MAIN_PROGRAM)_release: $(RELEASE_OBJ)
	$(LINK) $(RELEASE_LINKFLAGS) -o $@ $^ $(LINK_LIBS)

.PHONY: release
release: $(MAIN_PROGRAM)_release

.PHONY: debug
debug: $(MAIN_PROGRAM)

# how to compile each file
.SUFFIXES:
obj/%.o:
	$(CC) -c $(RELEASE_FLAGS) -o $@ $<
obj_debug/%.o:
	$(CC) -c $(DEBUG_FLAGS) -o $@ $<

# cleaning up
.PHONY: clean
clean:
	-rm -f obj/*.o $(MAIN_PROGRAM) obj_debug/*.o $(MAIN_PROGRAM)_release *core viewer

# dependencies are automatically generated
.PHONY: depend
depend:
	-mkdir obj
	-rm -f obj/depend
	$(foreach srcfile,$(SRC),$(DEPEND) -MM $(srcfile) -MT $(patsubst %.cpp,obj/%.o,$(notdir $(srcfile))) >> obj/depend;)
	-mkdir obj_debug
	-rm -f obj_debug/depend
	$(foreach srcfile,$(SRC),$(DEPEND) -MM $(srcfile) -MT $(patsubst %.cpp,obj_debug/%.o,$(notdir $(srcfile))) >> obj_debug/depend;)

-include obj/depend
-include obj_debug/depend

