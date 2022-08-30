// ---------------------------------------- //
// physics [C++ File]
// Written By: Thomas Bement
// Created On: 2022-08-29
// ---------------------------------------- //

// INCLUDES //
#include <iostream>

#include <string>

#include <vector>
using std::vector;

#include "planets.h"

// GLOBAL VARIABLES //
const float dt = 1;
const float G = 6.6743e-11; // [m.m.m/kg.s.s]

// FUNCTIONS //
void a_up(float dt, vector<planet> &planets) {
    for (int i = 0; i < planets.size(); i++) {
        for (int j = 0; j < i; j++) {
            vec3 displacement = planets[i].p.sub(planets[j].p);
            float distance = displacement.abs();
            vec3 G_inverse_square_ij = displacement.scale(G / (distance * distance * distance));
            planets[i].a = G_inverse_square_ij.scale(-planets[j].mass);
            planets[j].a = G_inverse_square_ij.scale(planets[i].mass);
        }
    }
    for (int i = 0; i < planets.size(); i++) {
        planets[i].v_up(dt);
        planets[i].p_up(dt);
    }
}

// MAIN //
int main() {
    std::vector<planet> planets = {};

    planet body_1(5.97219e24, vec3(0, 0, 0), vec3(0, 0, 0));
    planet body_2(7.34767309e22, vec3(3.844e8, 0, 0), vec3(0, 1022, 0));

    planets.push_back(body_1);
    planets.push_back(body_2);
    
    std::cout << planets[0].p << "   " << planets[1].p << std::endl;
    for (int i = 0; i < 1e9; i++) {
        a_up(dt, planets);
    }
    std::cout << planets[0].p << "   " << planets[1].p << std::endl;
    
    return 0;
}
