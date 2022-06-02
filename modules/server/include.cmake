if (APPLE OR WIN32)
  set(DEFAULT_MODULE ON)
else ()
  # Server is not available on Linux
  set(DEFAULT_MODULE OFF)
endif ()

set(OPENSPACE_DEPENDENCIES
  skybrowser
)
