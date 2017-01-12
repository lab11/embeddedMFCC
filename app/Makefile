# makefile for user application

TOCK_USERLAND_BASE_DIR = ../../tock/userland

C_SRCS += $(wildcard *.c)

CFLAGS += -D__CORTEX_M=0x04 -g 
CPPFLAGS += -DSTACK_SIZE=2048

# include userland master makefile. Contains rules and flags for actually
# 	building the application
include $(TOCK_USERLAND_BASE_DIR)/Makefile

.PHONY:
clean::
	rm -Rf $(BUILDDIR)
