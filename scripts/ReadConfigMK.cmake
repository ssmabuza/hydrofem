# Read config.mk
file(READ ${CMAKE_SOURCE_DIR}/config.mk CONFIG_MK_CONTENT)

# Extract variable definitions
string(REGEX MATCHALL "([A-Za-z0-9_]+) ?= ?(.*)\n" VAR_DEFS "${CONFIG_MK_CONTENT}")

# Loop through definitions and set variables
foreach(VAR_DEF ${VAR_DEFS})
    string(REGEX REPLACE "([A-Za-z0-9_]+) ?= ?(.*)\n" "\\1;\\2" VAR_DEF_CLEAN "${VAR_DEF}")
    list(GET VAR_DEF_CLEAN 0 VAR_NAME)
    list(GET VAR_DEF_CLEAN 1 VAR_VALUE)
    
    # Set the variable in CMake, avoiding overwriting existing variables
    if(NOT DEFINED ${VAR_NAME})
      set(${VAR_NAME} "${VAR_VALUE}" CACHE STRING "Imported from config.mk" FORCE)
    endif()
endforeach()

# Example usage: print the value of a variable
message(STATUS "Value of MY_VARIABLE: ${MY_VARIABLE}")
