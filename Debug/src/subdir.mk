################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/Atom.cpp \
../src/PhiDerivate.cpp \
../src/RadialDerivate.cpp \
../src/ThetaDerivate.cpp \
../src/trimero.cpp 

OBJS += \
./src/Atom.o \
./src/PhiDerivate.o \
./src/RadialDerivate.o \
./src/ThetaDerivate.o \
./src/trimero.o 

CPP_DEPS += \
./src/Atom.d \
./src/PhiDerivate.d \
./src/RadialDerivate.d \
./src/ThetaDerivate.d \
./src/trimero.d 


# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ -O0 -g3 -Wall -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


