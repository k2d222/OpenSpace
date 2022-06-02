if (APPLE OR WIN32)
  set(DEFAULT_MODULE ON)
else ()
  # Webgui is not available on Linux
  set(DEFAULT_MODULE OFF)
endif ()
