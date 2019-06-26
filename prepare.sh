#!/bin/bash
git clone https://github.com/alexkernphysiker/RectangularScintillator.git
git clone https://github.com/alexkernphysiker/math_h.git
cd math_h
	git checkout 6c507241e26e3ba435570234845aab75a8a0ad65
cd ..
git clone https://github.com/google/googletest.git
mkdir build
cd build
cmake ..
if make; then
	if math_h/math_h_test.exe --gtest_break_on_failure; then
		if RectangularScintillator/rectscin-test.exe --gtest_break_on_failure; then
			echo "Everythong is OK. You can run the calculations"
		fi
	fi
fi
cd ..
