#!/bin/bash
git clone https://github.com/alexkernphysiker/RectangularScintillator.git
git clone https://github.com/alexkernphysiker/math_h.git
cd math_h
#ToDo checkout
cd ..
git clone https://github.com/google/googletest.git
mkdir build
cd build
cmake ..
if make; then
math_h/math_h_test.exe
RectangularScintillator/rectscin-test.exe
fi
cd ..
