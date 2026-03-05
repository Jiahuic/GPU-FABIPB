CC ?= gcc

SRC_DIR ?= src
INC_DIR ?= include
BUILD_DIR ?= build
BIN ?= coulomb

# Prefer environment overrides so LAPACK/BLAS can come from the system,
# Homebrew, MKL, OpenBLAS, etc.
PROJECT_CPPFLAGS = -I$(INC_DIR)
CPPFLAGS ?=
CFLAGS ?= -g
LDFLAGS ?=
LAPACK_LIBS ?= -llapack
BLAS_LIBS ?= -lblas
MATH_LIBS ?= -lm

# Some platforms need libf2c explicitly for f2c-translated code.
# Example:
#   make F2C_LIBS=-lf2c
F2C_LIBS ?=

LIBS = $(LAPACK_LIBS) $(BLAS_LIBS) $(MATH_LIBS) $(F2C_LIBS)

HEADERS = $(INC_DIR)/gk.h $(INC_DIR)/gkGlobal.h $(INC_DIR)/f2c.h

SRCS = $(SRC_DIR)/coulomb.c \
       $(SRC_DIR)/gkSetup.c \
       $(SRC_DIR)/input.c \
       $(SRC_DIR)/gkGlobal.c \
       $(SRC_DIR)/expan.c \
       $(SRC_DIR)/moments.c \
       $(SRC_DIR)/numQuad.c \
       $(SRC_DIR)/fmm.c \
       $(SRC_DIR)/kernel.c \
       $(SRC_DIR)/precond_fmm.c \
       $(SRC_DIR)/gmres.c \
       $(SRC_DIR)/treecode.c

OBJS = $(SRCS:$(SRC_DIR)/%.c=$(BUILD_DIR)/%.o)

.PHONY: all clean print-config

all: $(BIN)

$(BUILD_DIR):
	mkdir -p $(BUILD_DIR)

$(BUILD_DIR)/%.o: $(SRC_DIR)/%.c $(HEADERS) Makefile | $(BUILD_DIR)
	$(CC) $(CPPFLAGS) $(PROJECT_CPPFLAGS) $(CFLAGS) -c -o $@ $<

$(BIN): $(OBJS)
	$(CC) $(LDFLAGS) -o $@ $(OBJS) $(LIBS)

print-config:
	@printf 'CC=%s\n' "$(CC)"
	@printf 'CPPFLAGS=%s\n' "$(CPPFLAGS)"
	@printf 'PROJECT_CPPFLAGS=%s\n' "$(PROJECT_CPPFLAGS)"
	@printf 'CFLAGS=%s\n' "$(CFLAGS)"
	@printf 'LDFLAGS=%s\n' "$(LDFLAGS)"
	@printf 'LAPACK_LIBS=%s\n' "$(LAPACK_LIBS)"
	@printf 'BLAS_LIBS=%s\n' "$(BLAS_LIBS)"
	@printf 'F2C_LIBS=%s\n' "$(F2C_LIBS)"

clean:
	rm -rf $(BUILD_DIR) $(BIN)
