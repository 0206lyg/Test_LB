# Application and private OpenLB core use the same official config/rules.
# Upstream config.mk, source files and existing core archive are not rewritten.
include $(OLB_ROOT)/config.mk
include $(OLB_ROOT)/rules.mk

.DEFAULT_GOAL := all
CORE_NAMES := communication/mpiManager communication/ompManager core/olbInit core/expr io/ostreamManager
CORE_OBJECTS := $(addprefix $(BUILD_DIR)/core/,$(addsuffix .o,$(CORE_NAMES)))
APP_OBJECT := $(BUILD_DIR)/slurry.o
INCLUDES := -I$(OLB_ROOT)/src -I$(OLB_ROOT)/external/zlib -I$(OLB_ROOT)/external/tinyxml2 -I$(BUILD_DIR)
ifeq ($(PARTICLE_SOLVER),petsc)
  INCLUDES += $(PETSC_CFLAGS)
  CXXFLAGS += -DSLURRY_USE_PETSC
  LDFLAGS += $(PETSC_LIBS)
endif

.PHONY: all
all: $(BUILD_DIR)/slurry

$(BUILD_DIR)/core/%.o: $(OLB_ROOT)/src/%.cpp
	@mkdir -p $(@D)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -MMD -MP -c $< -o $@

$(APP_OBJECT): $(BUILD_DIR)/slurry_application.cpp $(BUILD_DIR)/slurry_registry.h
	$(CXX) $(CXXFLAGS) $(INCLUDES) -MMD -MP -c $< -o $@

$(BUILD_DIR)/libolbcore.a: $(CORE_OBJECTS)
	$(AR) rc $@ $^

$(BUILD_DIR)/slurry: $(APP_OBJECT) $(BUILD_DIR)/libolbcore.a
	$(CXX) $(APP_OBJECT) -L$(BUILD_DIR) -lolbcore -L$(OLB_ROOT)/external/lib $(LDFLAGS) -o $@

-include $(CORE_OBJECTS:.o=.d) $(APP_OBJECT:.o=.d)
