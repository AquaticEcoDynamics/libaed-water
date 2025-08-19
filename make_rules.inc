###############################################################################
#                                                                             #
# make_rules.inc  - include rules for making libaed                           #
#                                                                             #
#  Developed by :                                                             #
#      AquaticEcoDynamics (AED) Group                                         #
#      School of Agriculture and Environment                                  #
#      The University of Western Australia                                    #
#                                                                             #
#      http://aquatic.science.uwa.edu.au/                                     #
#                                                                             #
#  Copyright 2013 - 2025 -  The University of Western Australia               #
#                                                                             #
#   AED is free software: you can redistribute it and/or modify               #
#   it under the terms of the GNU General Public License as published by      #
#   the Free Software Foundation, either version 3 of the License, or         #
#   (at your option) any later version.                                       #
#                                                                             #
#   AED is distributed in the hope that it will be useful,                    #
#   but WITHOUT ANY WARRANTY; without even the implied warranty of            #
#   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             #
#   GNU General Public License for more details.                              #
#                                                                             #
#   You should have received a copy of the GNU General Public License         #
#   along with this program.  If not, see <http://www.gnu.org/licenses/>.     #
#                                                                             #
###############################################################################


all: ${TARGET}

${libdir}/${OUTLIB}.${so_ext}: ${libdir}/${OUTLIB}.a
	$(F90) ${SHARED} -o $@.${SOVERS}.${VERS} ${OBJS} ${LDFLAGS} ${SOFLAGS}
	ln -sf ${OUTLIB}.${so_ext}.${SOVERS}.${VERS} $@
	ln -sf ${OUTLIB}.${so_ext}.${SOVERS}.${VERS} $@.${SOVERS}

${libdir}/$(OUTLIB).a: ${objdir} ${moddir} ${libdir} ${OBJS}
	ar rv $@ ${OBJS}
	ranlib $@

#${TARGET}: ${objdir} ${moddir} ${libdir} ${OBJS}
#	ar rv $@ ${OBJS}
#	ranlib $@

${objdir}/%.o: ${srcdir}/%.F90 ../libaed-water/${srcdir}/aed_core.F90 ../libaed-water/${incdir}/aed.h
	$(F90) $(FFLAGS) -c $< -o $@

${objdir}:
	@mkdir ${objdir}

${moddir}:
	@mkdir ${moddir}

${libdir}:
	@mkdir ${libdir}

clean: ${objdir}
	@touch ${objdir}/1.o 1.i90
	@/bin/rm ${objdir}/*.o *.i90

distclean: clean
	@/bin/rm -rf ${libdir} ${moddir} ${objdir}
