// ---------------------------------------- //
// physics [C++ File]
// Written By: Thomas Bement
// Created On: 2022-08-29
// ---------------------------------------- //

// INCLUDES //
#include <iostream>
#include <string>

// GLOBAL VARIABLES //


// FUNCTIONS //


// CLASSES //

class vec3 {
  public:
    float x;
    float y;
    float z;
    
    vec3(float xx, float yy, float zz)
      : x(xx), y(yy), z(zz) {}
    
    vec3 p(vec3 v) {
        return vec3(x+v.x, y+v.y, z+v.z);
    }
};


// MAIN //

int main() 
{
    return 0;
}
