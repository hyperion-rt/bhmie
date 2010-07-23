FC = ifort -module src/modules

OBJECTS =	src/types.f90 \
		src/material.f90 \
		src/distributions.f90 \
		src/bhmie.f90 \
		src/bhmie_caller.f90 \
		src/main.f90

bhmie:		$(OBJECTS)
		$(FC) $(OBJECTS) -o bin/bhmie
