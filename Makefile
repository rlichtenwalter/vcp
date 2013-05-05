##
# Copyright (C) 2013 by Ryan N. Lichtenwalter
# Email: rlichtenwalter@gmail.com
#
# This file is part of the Vertex Collocation Profiles code base.
#
# The Vertex Collocation Profiles code base is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#
# The Vertex Collocation Profiles code base is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with the Vertex Collocation Profiles code base. If not, see <http://www.gnu.org/licenses/>.
##


# USER PARAMETERS

# The maximum number of neighbors to allow a single vertex. This will probably be removed with no effect on interfaces or performance in a future release. For now, the only cost of setting it higher is a larger up-front allocation, but if it is too low, processing will fail.
MAX_NEIGHBORS := 16384

# BUILD SYSTEM AREA
SRCDIR := src
INCDIR := inc
ifeq ($(DEBUG),1)
	OBJDIR := debug_obj
	BINDIR := debug_bin
else
	OBJDIR := obj
	BINDIR := bin
endif
	   
BINS := vcp_generate vcp_map directed_to_undirected ell_2_pairs
SRCS := $(wildcard $(SRCDIR)/*.cpp)
HEADERS := $(wildcard $(INCDIR)/*.hpp)
OBJS := $(patsubst $(SRCDIR)/%,$(OBJDIR)/%,$(SRCS:.cpp=.o))
DEPS := $(OBJS:.o=.d)

VCP_INCLUDE := -I ./$(INCDIR)
BOOST_INCLUDE := -I ./lib/boost_1_53_0
TCLAP_INCLUDE := -I ./lib/tclap-1.2.1/include

CC := g++
COMMON_FLAGS := -Wall -Wextra -Werror -Wno-unused-local-typedefs -std=c++11 -pedantic $(VCP_INCLUDE) $(BOOST_INCLUDE) $(TCLAP_INCLUDE) -D MAX_NEIGHBORS=$(MAX_NEIGHBORS)
ifeq ($(DEBUG),1)
	CPP_FLAGS := $(COMMON_FLAGS) -Og -g
else
	CPP_FLAGS := $(COMMON_FLAGS) -O3 -D NDEBUG
endif

all: $(addprefix $(BINDIR)/,$(BINS)) $(BINDIR)

$(OBJDIR):
	@- mkdir -p $(OBJDIR) 2> /dev/null || true

$(BINDIR):
	@- mkdir -p $(BINDIR) 2> /dev/null || true

$(BINDIR)/%: $(OBJDIR)/%.o | $(BINDIR)
	$(CC) $(CPP_FLAGS) $< -o $@;

$(OBJDIR)/%.o: $(SRCDIR)/%.cpp | $(OBJDIR)
	$(CC) $(CPP_FLAGS) -MMD -MP -c $< -o $@;

-include $(DEPS)

.PHONY: all clean

clean:
	- rm -f $(OBJS);
	- rm -f $(DEPS);
	- rm -f $(addprefix $(BINDIR)/,$(BINS));
	- rmdir $(OBJDIR) $(BINDIR) 2> /dev/null || true
	
