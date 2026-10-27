CXX = g++
CPPFLAGS ?=
CXXFLAGS ?= -O2 -Wall -Wno-unused-result -Wno-unknown-pragmas
LDFLAGS ?=
LDLIBS ?=

SOURCE_DIR := source_for_Cxx_shared_library
TARGET := smoothing_on_sphere_Cxx_shared_library.so
SOURCES := $(SOURCE_DIR)/CC_smoothing_on_sphere_python_lib.cc \
           $(SOURCE_DIR)/CU_smoothing_on_sphere_code.cc \
           $(SOURCE_DIR)/CU_kdtree_with_index.cc

.PHONY: all clean

all: $(TARGET)

$(TARGET): $(SOURCES) Makefile
	@printf 'Compiling %s...\n' '$@'
	@printf '%s\n' '$(CXX) $(CPPFLAGS) $(CXXFLAGS) -std=c++11 -fopenmp -shared -fPIC $(LDFLAGS) -o $@ $(SOURCE_DIR)/CC_smoothing_on_sphere_python_lib.cc $(LDLIBS)'
	@if $(CXX) $(CPPFLAGS) $(CXXFLAGS) -std=c++11 -fopenmp -shared -fPIC $(LDFLAGS) -o $@ $(SOURCE_DIR)/CC_smoothing_on_sphere_python_lib.cc $(LDLIBS); then \
		printf 'Success: %s is ready to use.\n' '$@'; \
	else \
		status=$$?; \
		printf 'Error: compilation of %s failed. See the compiler messages above.\n' '$@' >&2; \
		exit $$status; \
	fi

clean:
	$(RM) $(TARGET)
