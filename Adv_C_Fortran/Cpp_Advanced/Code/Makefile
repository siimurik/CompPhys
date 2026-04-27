all: Driver ComplexNumberTestSuiteRunner LinearSystemTestSuiteRunner BvpOdeTestSuiteRunner 

### Single executable from Chapter 12:
Driver:	Driver.cpp BoundaryConditions.o BvpOde.o FiniteDifferenceGrid.o Vector.o Matrix.o LinearSystem.o
	g++ -O -o Driver $^

### The testing framework is a two-step process
# 1. Header to C++ main program via cxxtest generating script
#eg.
#ComplexNumberTestSuite.cpp: 	ComplexNumberTestSuite.hpp
#	cxxtestgen --have-eh --error-printer -o ComplexNumberTestSuite.cpp ComplexNumberTestSuite.hpp
%TestSuite.cpp: %TestSuite.hpp
	cxxtestgen --have-eh --error-printer -o $@ $<

# 2. C++ main program to executable - Then run the executable with -v "verbose trace"
ComplexNumberTestSuiteRunner:	ComplexNumberTestSuite.cpp ComplexNumber.o
	g++ -g -o ComplexNumberTestSuiteRunner $^\
	&& ./ComplexNumberTestSuiteRunner -v
LinearSystemTestSuiteRunner:	LinearSystemTestSuite.cpp Matrix.o Vector.o LinearSystem.o
	g++ -g -o LinearSystemTestSuiteRunner $^\
	&& ./LinearSystemTestSuiteRunner -v
BvpOdeTestSuiteRunner:	BvpOdeTestSuite.cpp BoundaryConditions.o BvpOde.o FiniteDifferenceGrid.o Vector.o Matrix.o LinearSystem.o
	g++ -g -o BvpOdeTestSuiteRunner $^\
	&& ./BvpOdeTestSuiteRunner -v

### Instructions for building the classes						
%.o:	%.cpp %.hpp
	g++ -g -c $<

clean:
	rm -f *.o *TestSuite.cpp
