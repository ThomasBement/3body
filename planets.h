// ---------------------------------------- //
// planets [C++ Header File]
// Written By: Thomas Bement
// Created On: 2022-08-30
// ---------------------------------------- //

// INCLUDES //
#include <iostream>
#include <cmath>

class vec3 {
  public:
    float x;
    float y;
    float z;
    
    vec3(float xx, float yy, float zz)
      : x(xx), y(yy), z(zz) {}
    
    vec3 plus(vec3 v) {
        return vec3(x+v.x, y+v.y, z+v.z);
    }

    vec3 sub(vec3 v) {
        return vec3(x-v.x, y-v.y, z-v.z);
    }

    vec3 scale(float a) {
        return vec3(a*x, a*y, a*z);
    }

    float abs() {
        return sqrt(x*x + y*y + z*z);
    }
};

std::ostream &operator<<(std::ostream &os, vec3 const &m) { 
    return os << "(" << m.x << "," << m.y << "," << m.z << ")";
}

class planet {       
    public:            
        float mass;
        vec3 p;
        vec3 v;
        vec3 a;
        planet(float mm, vec3 pp, vec3 vv) 
            : mass(mm), p(pp), v(vv), a(0.0, 0.0, 0.0) {}
    
    void p_up(float dt) {
        this->p = this->p.plus(this->v.scale(dt));
    }

    void v_up(float dt) {
        this->v = this->v.plus(this->a.scale(dt));
    }
};
