cmake_minimum_required(VERSION 3.20)
project(tinyrender)

set(CMAKE_CXX_STANDARD 20)

add_executable(tinyrender main.cpp
        softrender/TGA_Image.cpp
        softrender/TGA_Image.h
        softrender/ImageDraw.h

        softrender/Math/Vector.h
        softrender/Math/Quaternation.h
        softrender/Math/Quaternation.cpp
        softrender/Math/Matrix.h
        softrender/Math/Matrix.cpp
        softrender/Math/Color.h
        softrender/Math/Color.cpp
        softrender/Math/MathHelper.h
        softrender/Math/MathDefine.h
        softrender/Math/Math.h
        softrender/Math/Math.cpp

        softrender/Render/Model.cpp
        softrender/Render/Model.h)
